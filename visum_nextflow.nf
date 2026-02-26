#!/usr/bin/env nextflow
/*
 * ==========================================================
 *  viSUM (DSL2)
 *  Author: <Your Name>
 *  Project: <Project Name>
 *  Date: <Date>
 * ==========================================================
 */

nextflow.enable.dsl=2

#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process NORMALIZE_FASTA {

  tag "$prefix"

  // Per-sample outputs go here
  publishDir { "${params.outdir}/${prefix}/prep" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(fasta)

  output:
    tuple val(prefix), val(type), path("${prefix}.normalized.fasta"), path("${prefix}.header_map.tsv")

  script:
  """
  python3 - << 'PY'
  import sys

  prefix = ${prefix!r}
  in_fa  = ${fasta!r}
  out_fa = f"{prefix}.normalized.fasta"
  out_map= f"{prefix}.header_map.tsv"

  def fasta_reader(fp):
    header = None
    seq_chunks = []
    for line in fp:
      line = line.rstrip("\\n")
      if not line:
        continue
      if line.startswith(">"):
        if header is not None:
          yield header, "".join(seq_chunks)
        header = line[1:]  # keep full header (without >)
        seq_chunks = []
      else:
        seq_chunks.append(line.strip())
    if header is not None:
      yield header, "".join(seq_chunks)

  n = 0
  seen_new = set()

  with open(in_fa, "r", encoding="utf-8", errors="replace") as fin, \\
       open(out_fa, "w", encoding="utf-8") as fout, \\
       open(out_map, "w", encoding="utf-8") as fmap:

    # TSV header
    fmap.write("prefix\\tnew_id\\toriginal_id\\toriginal_header\\tlength\\n")

    for header, seq in fasta_reader(fin):
      if header is None:
        continue

      n += 1
      new_id = f"{prefix}__c{n:06d}"

      # Safety check (should never collide unless n resets)
      if new_id in seen_new:
        raise RuntimeError(f"Duplicate new_id generated: {new_id}")
      seen_new.add(new_id)

      # original_id: first token of header (common FASTA convention)
      original_id = header.split()[0] if header.strip() else ""

      # Write normalized fasta
      fout.write(f">{new_id}\\n")
      # wrap sequence to 80 chars for readability
      for i in range(0, len(seq), 80):
        fout.write(seq[i:i+80] + "\\n")

      # Mapping row
      fmap.write(f"{prefix}\\t{new_id}\\t{original_id}\\t{header}\\t{len(seq)}\\n")

  if n == 0:
    raise RuntimeError(f"No FASTA records found in {in_fa}")

  PY
  """
}

/*
  ORF prediction with pyrodigal-gv

  Input tuple:  (prefix, type, normalized_fasta, header_map)
  Output tuple: (prefix, type, normalized_fasta, header_map, proteins_faa)

  We also *save* the other outputs (genes + gff) for user review,
  but downstream we only require the .faa for later protein-based tools.
*/
process ORF_PREDICTION {

  tag "$prefix"

  conda 'bioconda::pyrodigal-gv'

  publishDir { "${params.outdir}/${prefix}/orf_prediction" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map)

  output:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path("${prefix}.protein.faa")
    path("${prefix}.orf.gff")
    path("${prefix}.gene.fna")

  script:
  """
  pyrodigal-gv -p meta -i ${norm_fasta} -a ${prefix}.protein.faa -d ${prefix}.gene.fna -f gff -o ${prefix}.orf.gff -q -j ${params.threads}
  """
}

process CHECKV_DB {

  tag "checkv_db"

  conda 'bioconda::checkv'

  publishDir { "${params.dbdir}/checkv" }, mode: 'copy'

  output:
    path("checkv_db")

  script:
  """
  set -euo pipefail
  mkdir -p checkv_db
  checkv download_database checkv_db
  """
}

process DEEP6_DB {

  tag "deep6_db"

  // git is needed to clone; python not strictly needed here
  conda 'conda-forge::git'

  publishDir { "${params.dbdir}/deep6" }, mode: 'copy'

  output:
    path("Deep6")

  script:
  """
  set -euo pipefail
  git clone --depth 1 https://github.com/janfelix/Deep6.git Deep6
  """
}

process GENOMAD_DB {

  tag "genomad_db"

  conda 'bioconda::genomad'

  // put the db in a stable location outside work/
  publishDir { "${params.dbdir}/genomad" }, mode: 'copy'

  output:
    path("genomad_db")

  script:
  """
  mkdir -p genomad_db
  genomad download-database genomad_db
  """
}

process RUN_CHECKV {

  tag "$prefix"

  conda 'bioconda::checkv'

  publishDir { "${params.outdir}/${prefix}/checkv" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path checkv_db

  output:
    tuple val(prefix), val("checkv"), path("${prefix}.checkv_quality_summary.tsv")

  script:
  """
  set -euo pipefail

  outdir="${prefix}.checkv.review"
  mkdir -p "\$outdir"

  checkv end_to_end ${norm_fasta} "\$outdir" -t ${params.threads} -d ${checkv_db}

  # If checkv produced the summary, copy it to stable name; else create header-only
  if [ -f "\$outdir/quality_summary.tsv" ]; then
    cp "\$outdir/quality_summary.tsv" ${prefix}.checkv_quality_summary.tsv
  else
    echo -e "contig_id\\tcheckv_quality\\tcompleteness\\tcontamination\\tmethod\\twarnings" > ${prefix}.checkv_quality_summary.tsv
  fi
  """
}

process RUN_DEEP6 {

  tag "$prefix"

  /*
    Deep6 README lists deps like numpy/pandas/h5py/biopython/scipy/keras/tensorflow/scikit-learn. :contentReference[oaicite:3]{index=3}
    You may want to pin versions later once you test on your system.
  */
  conda 'conda-forge::python=3.10 conda-forge::numpy conda-forge::pandas conda-forge::h5py conda-forge::biopython conda-forge::scipy conda-forge::scikit-learn conda-forge::tensorflow conda-forge::keras'

  publishDir { "${params.outdir}/${prefix}/deep6" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path deep6_db

  output:
    tuple val(prefix), val("deep6"), path("${prefix}.${params.deep6_minlen}bp_deep6_evidence.csv")

  when:
    type == 'rna'

  script:
  """
  set -euo pipefail

  # -------- model dir (user override or default to repo Models/) ----------
  MODEL_DIR="${params.deep6_model}"
  if [ -z "${MODEL_DIR}" ] || [ "${MODEL_DIR}" = "null" ]; then
    MODEL_DIR="${deep6_db}/Models"
  fi

  outdir="${prefix}.deep6.review"
  mkdir -p "$outdir"

  # -------- run deep6 ----------
  python3 ${deep6_db}/Master/deep6.py -i ${norm_fasta} -l ${params.deep6_minlen} -m "$MODEL_DIR" -o "$outdir"

# -------- locate prediction output ----------
  pred_file=$(ls "$outdir"/*_predict_${params.deep6_minlen}bp_deep6.txt 2>/dev/null || true)

  RAW_TSV="${prefix}.${params.deep6_minlen}bp_deep6_scores_raw.tsv"
  RAW_CSV="${prefix}.${params.deep6_minlen}bp_deep6_scores_raw.csv"
  CLEAN_CSV="${prefix}.${params.deep6_minlen}bp_deep6_scores_clean.csv"
  EVID="${prefix}.${params.deep6_minlen}bp_deep6_evidence.csv"
  CONF="${prefix}.${params.deep6_minlen}bp_deep6_confident_predictions.csv"

  # Always create evidence header (so Nextflow output exists)
  echo "seqid,d__Domain,r__Realm,k__Kingdom,p__Phylum,c__Class,o__Order,f__Family,g__Genus,s__Species" > "$EVID"

  # If deep6 produced predictions, copy them; else create header-only raw.tsv
  if [ -n "$pred_file" ] && [ -s "$pred_file" ]; then
    cp "$pred_file" "$RAW_TSV"
  else
    # Deep6 standard header per repo output: name length duplo euk mono pro ribo vari
    printf "name\tlength\tduplo\teuk\tmono\tpro\tribo\tvari\n" > "$RAW_TSV"
  fi

  # Convert tsv -> csv (overwrite, don't append)
  sed 's/\t/,/g' "$RAW_TSV" > "$RAW_CSV"

  # If RAW_CSV has more than header, run your scoring/cleaning
  # (NR>1 means at least 1 data row)
  if awk 'NR>1{exit 0} END{exit 1}' "$RAW_CSV"; then
    python3 ${projectDir}/bin/deep6_compare_all_vs_allv0.1.py "$RAW_CSV" "${prefix}_deep6_sanity_check.csv" "$CLEAN_CSV"
    # Build confident predictions + evidence
    # Assumes CLEAN_CSV format: seqid,length,score,realm,flag (as you described)
    echo "seqid,length,score,realm,flag" > "$CONF"
    # pull unique seqids (skip header)
    awk -F "," 'NR>1{print $1}' "$CLEAN_CSV" | sort -u > sample.list
    while read -r x; do
      top=$(grep -w "$x" "$CLEAN_CSV" | awk -F "," '{print $3}' | head -n 1)
      # guard if missing
      if [ -z "$top" ]; then
        continue
      fi
      # if top > 0.7 keep it
      if (( $(echo "$top > 0.7" | bc -l) )); then
        grep -w "$x" "$CLEAN_CSV" >> "$CONF"
        # realm is column 4
        realm=$(grep -w "$x" "$CLEAN_CSV" | awk -F "," '{print $4}' | head -n 1)
        # Map realm -> taxonomy line using your mapping file
        # deep6lines.txt should include entries like: ribo;d__Viruses,r__Riboviria,k__unclassified,...
        line=$(grep -w "^${realm};" ${projectDir}/bin/deep6lines.txt | awk -F ";" '{print $2}' | head -n 1)
        # if mapping missing, default to viruses/unclassified
        if [ -z "$line" ]; then
          line="d__Viruses,r__unclassified,k__unclassified,p__unclassified,c__unclassified,o__unclassified,f__unclassified,g__unclassified,s__unclassified"
        fi
        echo "${x},${line}" >> "$EVID"
      fi
    done < sample.list
  else
    # no data rows -> evidence stays header-only
    :
  fi
  """
}

process RUN_GENOMAD {

  tag "$prefix"

  conda 'bioconda::genomad'

  publishDir { "${params.outdir}/${prefix}/genomad" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path genomad_db

  output:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path("${prefix}.genomad.review"),
          path("${prefix}.genomad_virus_summary.tsv"),
          path("${prefix}.genomad_plasmid_summary.tsv")

  script:
  """
  set -euo pipefail
  base="${norm_fasta.baseName}"
  outdir="${prefix}.genomad.review"
  genomad end-to-end --cleanup --threads ${params.threads} ${norm_fasta} "\$outdir" ${genomad_db}
  summary_dir="\$outdir/\${base}_summary"
  virus_src="\$summary_dir/\${base}_virus_summary.tsv"
  plasmid_src="\$summary_dir/\${base}_plasmid_summary.tsv"
  # Fail early with a helpful message if outputs aren't where we expect.
  test -s "\$virus_src"   || { echo "[ERROR] Missing virus summary: \$virus_src" >&2; ls -R "\$outdir" >&2; exit 1; }
  test -s "\$plasmid_src" || { echo "[ERROR] Missing plasmid summary: \$plasmid_src" >&2; ls -R "\$outdir" >&2; exit 1; }
  # Copy into stable filenames so Nextflow outputs are simple/robust
  cp "\$virus_src"   ${prefix}.genomad_virus_summary.tsv
  cp "\$plasmid_src" ${prefix}.genomad_plasmid_summary.tsv
  """
}

process GENOMAD_PROCESSING {

  tag "$prefix"

  conda 'bioconda::python=3.11 pandas'

  publishDir { "${params.outdir}/${prefix}/genomad" }, mode: 'copy'

  input:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path(genomad_review),
          path(virus_summary),
          path(plasmid_summary)
    path ictv_csv

  output:
    tuple val(prefix), val("genomad"), path("${prefix}.genomad.vsum.csv")

  script:
  """
  set -euo pipefail

  python3 ${projectDir}/bin/fixgenomadv0.2.py --genomad ${virus_summary} --ictv ${params.ictv_csv} --out ${prefix}.genomad.vsum.csv --fallback ictv

  if [ ! -s ${prefix}.genomad.vsum.csv ]; then
    echo "seqid,d__Domain,r__Realm,k__Kingdom,p__Phylum,c__Class,o__Order,f__Family,g__Genus,s__Species" > ${prefix}.genomad.vsum.csv
  fi
  """
}


workflow {

  /*
    -----------------------------
    Decide mode + validate inputs
    -----------------------------
  */
  def singleMode = (params.input != null)
  def multiMode  = (params.prefix_many != null)

  if( singleMode && multiMode )
    error "Choose ONE mode: (1) --input --prefix --type  OR  (2) --prefix_many --indir"

  if( !singleMode && !multiMode )
    error "Provide inputs using: (1) --input --prefix --type  OR  (2) --prefix_many --indir"

  Channel ch_samples

  /*
    -----------------------------
    SINGLE MODE
    emits one tuple: (prefix, type, fasta)
    -----------------------------
  */
  if( singleMode ) {

    if( !params.prefix )
      error "Single mode requires --prefix"
    if( !params.type )
      error "Single mode requires --type (dna|rna)"

    def t = params.type.toString().trim().toLowerCase()
    if( !(t in ['dna','rna']) )
      error "Invalid --type '${params.type}'. Must be dna or rna."

    def fasta = file(params.input)
    if( !fasta.exists() )
      error "Input FASTA not found: ${params.input}"

    ch_samples = Channel.of( tuple(params.prefix, t, fasta) )
  }

  /*
    -----------------------------
    MULTI MODE
    reads CSV + indir
    emits many tuples: (prefix, type, fasta)
    -----------------------------
  */
  if( multiMode ) {

    if( !params.indir )
      error "Multi mode requires --indir <directory containing FASTA files>"

    def baseDir = file(params.indir)
    if( !baseDir.exists() || !baseDir.isDirectory() )
      error "--indir must be an existing directory: ${params.indir}"

    def mapfile = file(params.prefix_many)
    if( !mapfile.exists() )
      error "prefix_many CSV not found: ${params.prefix_many}"

    ch_samples = Channel
      .fromPath(mapfile)
      .splitCsv(header: true, sep: ',')   // expects header: prefix,type,fasta
      .map { row ->

        def prefix = row.prefix?.toString()?.trim()
        def type   = row.type?.toString()?.trim()?.toLowerCase()
        def rel    = row.fasta?.toString()?.trim()

        if( !prefix || !type || !rel )
          error "CSV must have columns prefix,type,fasta with non-empty values. Bad row: ${row}"

        if( !(type in ['dna','rna']) )
          error "Invalid type for prefix=${prefix}: '${type}' (must be dna or rna)"

        def fasta = file("${baseDir}/${rel}")
        if( !fasta.exists() )
          error "FASTA not found for prefix=${prefix}: ${baseDir}/${rel}"

        tuple(prefix, type, fasta)
      }
  }

  /*
    -----------------------------
    Debug: print what will run
    -----------------------------
  */
  ch_samples.view { p, t, f -> "SAMPLE=${p}\tTYPE=${t}\tFASTA=${f}" }

// CHECKV DB
if( !file(params.ictv_csv).exists() )
    error "ICTV CSV not found: ${params.ictv_csv}"

def CHECKV_DB_PATH = file("${params.dbdir}/checkv/checkv_db")

Channel ch_checkv_db

if( CHECKV_DB_PATH.exists() ) {
  ch_checkv_db = Channel.value(CHECKV_DB_PATH)
} else {
  ch_checkv_db = CHECKV_DB().map { it -> CHECKV_DB_PATH }
}

//DEEP6 DB / SETUP
def DEEP6_DB_PATH = file(params.deep6_dir)

Channel ch_deep6_db

if( DEEP6_DB_PATH.exists() ) {
  ch_deep6_db = Channel.value(DEEP6_DB_PATH)
} else {
  ch_deep6_db = DEEP6_DB().map { it -> DEEP6_DB_PATH }
}

// GENOMAD DB
def GENOMAD_DB_PATH = file("${params.dbdir}/genomad/genomad_db")

Channel ch_genomad_db

if( GENOMAD_DB_PATH.exists() ) {
  // DB already present: just point to it
  ch_genomad_db = Channel.value(GENOMAD_DB_PATH)
} else {
  // DB missing: run download once
  ch_genomad_db = GENOMAD_DB().map { it -> GENOMAD_DB_PATH }
}

  // ---------- NORMALIZE ----------
  ch_norm = NORMALIZE_FASTA(ch_samples)

  // Debug: show outputs of normalization
  ch_norm.view { p, t, nf, map ->
    "NORM   sample=${p}\ttype=${t}\tnorm_fasta=${nf.name}\tmap=${map.name}"
  }

// ---------- ORF PREDICTION ----------
(ch_orf, ch_orf_gff, ch_orf_fna) = ORF_PREDICTION(ch_norm)

// Debug: show the *main* tuple stream that will feed downstream tools
ch_orf.view { p, t, nf, map, faa ->
  "ORF    sample=${p}\ttype=${t}\tfaa=${faa.name}\tnorm_fasta=${nf.name}"
}

def branched = ch_orf.branch(
  dna: { p, t, nf, hm, faa -> t == 'dna' },
  rna: { p, t, nf, hm, faa -> t == 'rna' }
)

// CHECKV PROCESSES
ch_checkv_ev = RUN_CHECKV(ch_orf, ch_checkv_db)

// DEEP6 PROCESSES


// GENOMAD PROCESSES
ch_genomad_raw = RUN_GENOMAD(ch_orf, ch_genomad_db)
ch_genomad_raw.view { p,t,nf,hm,faa,review,vs,ps ->
  "GENOMAD sample=${p}\tvirus_summary=${vs.name}\tplasmid_summary=${ps.name}"
}
ch_genomad_ev = GENOMAD_PROCESSING(ch_genomad_raw, ictv_csv_ch)

// BELOW IS THE TEMPLATE TO MERGE EVIDENCE FROM ALL PROGRAMS RAN ON FASTA FILE

// core tuple from ORF stage: (prefix, type, norm_fasta, header_map, proteins_faa)
ch_core = ch_orf.map { p, t, nf, hm, faa -> tuple(p, t, nf, hm, faa) }

// gather all evidence records (prefix, tool, file)
ch_all_evidence = ch_genomad_ev
  // .mix(ch_checkv_ev)
  // .mix(ch_virsorter_ev)
  // .mix(ch_diamond_ev)
  .groupTuple(by: 0)

ch_all_evidence.view { prefix, tools, files ->
  "EVIDENCE prefix=${prefix}\ttools=${tools}\tfiles=${files.collect{ it.name }}"
}

// convert grouped evidence lists -> map
ch_ev_map = ch_all_evidence.map { prefix, tools, files ->
  def m = [:]
  tools.eachWithIndex { t, i -> m[t] = files[i] }
  tuple(prefix, m)
}

ch_ev_map.view { prefix, m ->
  "EV_MAP  prefix=${prefix}\tkeys=${m.keySet().sort()}"
}

// join core + evidence map -> viSUM inputs
ch_vsum_in = ch_core.join(ch_ev_map)
  .map { core, ev ->
    def (p, t, nf, hm, faa) = core
    def (_, m) = ev
    tuple(p, t, nf, hm, faa, m)
  }

}
