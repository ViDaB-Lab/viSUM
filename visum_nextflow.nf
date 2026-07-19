#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
 * ==========================================================
 *  viSUM (DSL2)
 *  Author: <Your Name>
 *  Project: <Project Name>
 *  Date: <Date>
 * ==========================================================
 */




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

  conda 'bioconda::checkv conda-forge::rsync'

  output:
    path("checkv_db")  // symlink to final DB dir

  script:
  """
  set -euo pipefail

  DB_DEST="${params.dbdir}/checkv/checkv_db"
  SENTINEL="\${DB_DEST}/.db_complete"

  mkdir -p "\${DB_DEST}"

  # Already installed
  if [[ -f "\${SENTINEL}" ]]; then
    ln -sfn "\${DB_DEST}" checkv_db
    exit 0
  fi

  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  # Download into temp first to avoid partial DB in final location
  checkv download_database "\$tmpdir/checkv_db"

  # Sync into place (robust across filesystems)
  rsync -a --delete "\$tmpdir/checkv_db/" "\${DB_DEST}/"

  # Mark complete
  date -Iseconds > "\${SENTINEL}"

  # Emit symlink output
  ln -sfn "\${DB_DEST}" checkv_db
  """
}

process DEEP6_DB {

  tag "deep6_db"

  conda 'conda-forge::git conda-forge::rsync'

  output:
    path("Deep6")   // symlink to final Deep6 repo dir

  script:
  """
  set -euo pipefail

  DB_DEST="${params.deep6_dir}"
  SENTINEL="\${DB_DEST}/.db_complete"

  mkdir -p "\${DB_DEST}"

  # Already installed
  if [[ -f "\${SENTINEL}" ]]; then
    ln -sfn "\${DB_DEST}" Deep6
    exit 0
  fi

  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  # Clone into temp first (avoid partial repo in final location)
  git clone --depth 1 https://github.com/janfelix/Deep6.git "\$tmpdir/Deep6"

  # Sync into place
  rsync -a --delete "\$tmpdir/Deep6/" "\${DB_DEST}/"

  # Mark complete
  date -Iseconds > "\${SENTINEL}"

  # Emit symlink output
  ln -sfn "\${DB_DEST}" Deep6
  """
}

process GENOMAD_DB {

  tag "genomad_db"

  conda 'bioconda::genomad conda-forge::rsync'

  output:
    path("genomad_db")  // symlink to final DB dir

  script:
  """
  set -euo pipefail

  DB_DEST="${params.dbdir}/genomad/genomad_db"
  SENTINEL="\${DB_DEST}/.db_complete"

  mkdir -p "\${DB_DEST}"

  # If already complete, just emit the symlink output
  if [[ -f "\${SENTINEL}" ]]; then
    ln -sfn "\${DB_DEST}" genomad_db
    exit 0
  fi

  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  # Download into temp first to avoid partial DB in final location
  genomad download-database "\$tmpdir/genomad_db"

  # Move into place (robust across filesystems)
  rsync -a --delete "\$tmpdir/genomad_db/" "\${DB_DEST}/"

  # Mark complete only after successful sync
  date -Iseconds > "\${SENTINEL}"

  # Emit symlink output for Nextflow
  ln -sfn "\${DB_DEST}" genomad_db
  """
}

process VIRSORTER2_DB {

  tag "virsorter2_db"

  conda 'bioconda::virsorter=2.2.4 conda-forge::rsync'

  output:
    path("db")  // symlink to real db dir

  script:
  """
  set -euo pipefail

  DB_DEST="${params.vs2_dir}"
  SENTINEL="\${DB_DEST}/.db_complete"

  mkdir -p "\${DB_DEST}"

  if [[ -f "\${SENTINEL}" ]]; then
    ln -sfn "\${DB_DEST}" db
    exit 0
  fi

  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  virsorter setup -d "\$tmpdir/db" -j ${params.threads}

  rsync -a --delete "\$tmpdir/db/" "\${DB_DEST}/"
  date -Iseconds > "\${SENTINEL}"

  ln -sfn "\${DB_DEST}" db
  """
}

process CENOTETAKER3_DB {

  tag "cenote_taker3_db"

  conda 'bioconda::cenote-taker3=3.* conda-forge::rsync'

  output:
    path("db")   // a symlink pointing to the real DB directory

  script:
    def hhFlags = params.ct3_use_hhsuite
      ? "--hhCDD T --hhPFAM T --hhPDB T"
      : ""

  """
  set -euo pipefail

  DB_DEST="${params.ct3_dir}"
  SENTINEL="\${DB_DEST}/.db_complete"

  mkdir -p "\${DB_DEST}"

  # If already complete, just emit the symlink output
  if [[ -f "\${SENTINEL}" ]]; then
    ln -sfn "\${DB_DEST}" db
    exit 0
  fi

  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  # Download into a temp dir first, then move into place to avoid partial installs
  get_ct3_dbs -o "\$tmpdir/db" \
    --hmm T \
    --hallmark_tax T \
    --refseq_tax T \
    --mmseqs_cdd T \
    --domain_list T \
    ${hhFlags}

  # Move into place (rsync is safer than mv across filesystems)
  rsync -a --delete "\$tmpdir/db/" "\${DB_DEST}/"

  date -Iseconds > "\${SENTINEL}"

  ln -sfn "\${DB_DEST}" db
  """
}

process VITAP_DB {

  tag "vitap_db:${params.vitap_db_label}"

  // VITAP provides its own CLI; include rsync for safe atomic copy
  conda 'bioconda::vitap=1.10 conda-forge::rsync'

  output:
    path("DB_${params.vitap_db_label}")

  script:
  """
  set -euo pipefail

  LABEL="${params.vitap_db_label}"
  VMR="${params.vitap_vmr_csv}"

  DBROOT="${params.vitap_dir}"
  DEST="\${DBROOT}/DB_\${LABEL}"
  SENTINEL="\${DEST}/.db_complete"

  mkdir -p "\${DBROOT}"

  # If already complete, just emit a symlink as the process output
  if [[ -f "\${SENTINEL}" ]]; then
    ln -s "\${DEST}" "DB_\${LABEL}"
    exit 0
  fi

  # Build in a temp directory first to avoid leaving partial DBs in the shared dbdir
  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  cp "\${VMR}" "\$tmpdir/vmr.csv"
  cd "\$tmpdir"

  # Non-interactive: auto-continue at the Y/N prompt
  printf "Y\\n" | VITAP upd --vmr vmr.csv -o VMR_reformat.csv -d "\${LABEL}"

  # VITAP should produce DB_<LABEL> in the working dir
  [[ -d "DB_\${LABEL}" ]] || { echo "[ERROR] Expected DB_\${LABEL} was not created." >&2; ls -lah >&2; exit 1; }

  # Copy into shared dbdir atomically/safely
  mkdir -p "\${DEST}"
  rsync -a --delete "DB_\${LABEL}/" "\${DEST}/"

  date -Iseconds > "\${SENTINEL}"

  # Emit a stable path as the process output
  ln -s "\${DEST}" "DB_\${LABEL}"
  """
}

process GIANTHUNTER_DB {

  tag "gianthunter_db"

  conda 'conda-forge::wget conda-forge::unzip conda-forge::rsync'

  output:
    path("gianthunter_db")

  script:
  """
  set -euo pipefail

  DB_DEST="${params.gianthunter_dir}"
  SENTINEL="\${DB_DEST}/.db_complete"

  mkdir -p "\${DB_DEST}"

  # If DB is already complete, just expose it to Nextflow
  if [[ -f "\${SENTINEL}" ]]; then
    ln -sfn "\${DB_DEST}" gianthunter_db
    exit 0
  fi

  tmpdir=\$(mktemp -d)
  trap 'rm -rf "\$tmpdir"' EXIT

  cd "\$tmpdir"
  wget -O gianthunter_db_v1.zip https://github.com/FuchuanQu/GiantHunter/releases/download/v2.0/gianthunter_db_v1.zip
  unzip gianthunter_db_v1.zip

  # Sync extracted DB into final destination
  rsync -a --delete "\$tmpdir/gianthunter_db_v1/" "\${DB_DEST}/"

  # Mark completion only after successful sync
  date -Iseconds > "\${SENTINEL}"

  # Expose stable path to Nextflow
  ln -sfn "\${DB_DEST}" gianthunter_db
  """
}

process VICAT_DB {

  tag "vicat_db:${params.vicat_db_label}"

  // VITAP provides its own CLI; include rsync for safe atomic copy
  conda 'bioconda::vitap=1.10 conda-forge::rsync'

  output:
    path("DB_${params.vicat_db_label}")

  script:
  """
  set -euo pipefail

  LABEL="${params.vicat_db_label}"

  DBROOT="${params.vicat_dir}"
  DEST="\${DBROOT}/DB_\${LABEL}"
  SENTINEL="\${DEST}/.db_complete"

  mkdir -p "\${DBROOT}"

  # If already complete, just emit a symlink as the process output
  if [[ -f "\${SENTINEL}" ]]; then
    ln -s "\${DEST}" "DB_\${LABEL}"
    exit 0
  fi

  # Build in a temp directory first to avoid leaving partial DBs in the shared dbdir
  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  cd "\$tmpdir"

  #download the proteins fasta file
  wget -c https://www.meta-virome.org/Data/Downloads/IMGVR5_UViG.faa.gz
  pigz -p ${params.threads} -d IMGVR5_UViG.faa.gz
  
  #download metadata file
  wget -c https://www.meta-virome.org/DownloadUvigMetadata 
  mv DownloadUvigMetadata ./DownloadUvigMetadata.gz
  gunzip DownloadUvigMetadata.gz

  # Copy into shared dbdir atomically/safely
  mkdir -p "\${DEST}"
  rsync -a --delete "DB_\${LABEL}/" "\${DEST}/"

  date -Iseconds > "\${SENTINEL}"

  # Emit a stable path as the process output
  ln -s "\${DEST}" "DB_\${LABEL}"
  """
}

process VCONTACT3_DB {

  tag "vcontact3_db"

  conda 'bioconda::vcontact3 conda-forge::rsync'

  output:
    path("vcontact3_db")

  script:
  """
  set -euo pipefail

  DB_DEST="${params.vcontact3_db_dir}"
  SENTINEL="\${DB_DEST}/.db_complete"

  mkdir -p "\${DB_DEST}"

  if [[ -f "\${SENTINEL}" ]]; then
    ln -sfn "\${DB_DEST}" vcontact3_db
    exit 0
  fi

  tmpdir="\$(mktemp -d)"
  trap 'rm -rf "\$tmpdir"' EXIT

  vcontact3 prepare_databases --get-version latest --set-location "\$tmpdir/db"

  rsync -a --delete "\$tmpdir/db/" "\${DB_DEST}/"

  if [[ ! -d "\${DB_DEST}" ]]; then
    echo "ERROR: vConTACT3 database setup failed" >&2
    exit 1
  fi

  date -Iseconds > "\${SENTINEL}"

  ln -sfn "\${DB_DEST}" vcontact3_db
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

  MODEL_DIR="${params.deep6_model}"
  if [ -z "${MODEL_DIR}" ] || [ "${MODEL_DIR}" = "null" ]; then
    MODEL_DIR="${deep6_db}/Models"
  fi
  outdir="${prefix}.deep6.review"
  mkdir -p "$outdir"

  #run deep6 
  python3 ${deep6_db}/Master/deep6.py -i ${norm_fasta} -l ${params.deep6_minlen} -m "$MODEL_DIR" -o "$outdir"

  # locate prediction output 
  pred_file=$(ls "$outdir"/*_predict_${params.deep6_minlen}bp_deep6.txt 2>/dev/null || true)

  # setting file names
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

process RUN_VIRSORTER2 {

  tag "$prefix"

  conda 'bioconda::virsorter=2.2.4'

  publishDir { "${params.outdir}/${prefix}/virsorter2" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path vs2_db

  output:
    tuple val(prefix), val("virsorter2"), path("${prefix}.virsorter2_evidence.csv")

  script:
  """
  set -euo pipefail

  outdir="${prefix}.virsorter2.review"
  mkdir -p "\$outdir"

  # Choose groups based on type (RNA-only vs broader DNA set)
  if [ "${type}" = "rna" ]; then
    GROUPS="${params.vs2_groups_rna}"
  else
    GROUPS="${params.vs2_groups_dna}"
  fi

  # VirSorter2 quick run pattern (run workflow, writes final-viral-score.tsv) 
  virsorter run -w "\$outdir" -i ${norm_fasta} --min-length ${params.vs2_min_length} --min-score ${params.vs2_min_score} --include-groups "\$GROUPS" -j ${params.threads} all

  # Copy the key output to a stable name for evidence aggregation :contentReference[oaicite:4]{index=4}
  if [ -f "\$outdir/final-viral-score.tsv" ]; then
    cp "\$outdir/final-viral-score.tsv" ${prefix}.virsorter2.final-viral-score.tsv
  else
    # Header-only safety net if VS2 produced nothing (keeps pipeline unbreakable)
    printf "seqname\\tmax_score\\tmax_group\\n" > ${prefix}.virsorter2.final-viral-score.tsv
  fi

  RAW="${prefix}.virsorter2.final-viral-score.tsv"
  EVID="${prefix}.virsorter2_evidence.csv"

  echo "seqid,vs2_tag,vs2_score,vs2_group,vs2_is_viral" > "\$EVID"

  if [ -s "\$RAW" ]; then
    awk -F "\t" '
    NR>1 {
      split(\$1, a, /\|\|/)
      seqid = a[1]
      tag   = (length(a)>1 ? a[2] : "full")
      score = \$7
      group = \$8
      print seqid","tag","score","group",1
  }' "\$RAW" >> "\$EVID"
  fi
  """
}

process RUN_CENOTETAKER3 {

  tag "$prefix"

  conda 'bioconda::cenote-taker3=3.*'

  publishDir { "${params.outdir}/${prefix}/ct3" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path ct3_db

  output:
    tuple val(prefix),
          val(type),
          path("${prefix}.ct3.review"),
          path("${prefix}.ct3_virus_summary.tsv")

  script:
  """
  set -euo pipefail

  # Point CT3 to the shared database directory
  export CENOTE_DBS="${ct3_db}"

  outdir="${prefix}.ct3.review"
  mkdir -p "\$outdir"
  run_title="${prefix}"

  cenote-taker3 -c ${norm_fasta} -r "\$run_title" -p ${params.prophage} -t ${params.threads} --caller ${params.ct3_caller} -hh ${params.ct3_hh} --taxdb ${params.ct3_taxdb} -db ${params.ct3_domaindb} --minimum_length_circular ${params.ct3_minlen_circ} --circ_minimum_hallmark_genes ${params.ct3_cir_minhall} --minimum_length_linear ${params.ct3_minlen_lin} --lin_minimum_hallmark_genes ${params.ct3_lin_minhall} -o "\$outdir" --genbank ${params.ct3_genb}

  # Locate virus summary
  vs="\$outdir/\${run_title}_virus_summary.tsv"
  test -f "\$vs" || { echo "[ERROR] Missing CT3 virus summary: \$vs" >&2; ls -R "\$outdir" >&2; exit 1; }

  cp "\$vs" ${prefix}.ct3_virus_summary.tsv
  """
}

process PROCESS_CENOTETAKER3 {

  tag "$prefix"

  conda 'bioconda::python=3.11 pandas'

  publishDir { "${params.outdir}/${prefix}/ct3" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(ct3_review), path(ct3_virus_summary)

  output:
    tuple val(prefix), val("ct3"), path("${prefix}.ct3_evidence.csv")

  script:
  """
  set -euo pipefail

  python3 ${projectDir}/bin/fixct3v0.2.py --ct3 ${ct3_virus_summary} --ictv ${params.ictv_csv} --out ${prefix}.ct3_evidence.csv --fallback ictv

  # Safety net
  if [ ! -s ${prefix}.ct3_evidence.csv ]; then
    echo "seqid,d__Domain,r__Realm,k__Kingdom,p__Phylum,c__Class,o__Order,f__Family,g__Genus,s__Species,ct3_virion_hallmark_count,ct3_rep_hallmark_count,ct3_RDRP_hallmark_count" > ${prefix}.ct3_evidence.csv
  fi
  """
}

process RUN_VITAP {

  tag "$prefix"

  conda 'bioconda::vitap=1.10 conda-forge::rsync'

  publishDir { "${params.outdir}/${prefix}/vitap" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path vitap_db

  output:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path("${prefix}.vitap.review"),
          path("${prefix}.vitap_best_determined_lineages.tsv")

  script:
  """
  set -euo pipefail

  outdir="${prefix}.vitap.review"
  mkdir -p "\$outdir"

  # Run VITAP assignment
  VITAP assignment -i ${norm_fasta} -d ${vitap_db} -o "\$outdir"

  # VITAP writes: <outdir>/best_determined_lineages.tsv
  raw="\$outdir/best_determined_lineages.tsv"

  # Always produce the file (even if empty) so downstream doesn't break
  if [[ -s "\$raw" ]]; then
    cp "\$raw" ${prefix}.vitap_best_determined_lineages.tsv
  else
    # create empty-but-valid tsv with header
    echo -e "Genome_ID\\tlineage\\tlineage_score/participation_index\\tConfidence_level" > ${prefix}.vitap_best_determined_lineages.tsv
  fi
  """
}

process PROCESS_VITAP {

  tag "$prefix"

  conda 'bioconda::python=3.11 pandas'

  publishDir { "${params.outdir}/${prefix}/vitap" }, mode: 'copy'

  input:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path(vitap_review),
          path(vitap_lineages)
    path ictv_csv

  output:
    tuple val(prefix),
          val("vitap"),
          path("${prefix}.vitap_evidence.csv")

  script:
  """
  set -euo pipefail

  OUT="${prefix}.vitap_evidence.csv"

  # Always create output with header so pipeline never breaks
  echo "seqid,d__Domain,r__Realm,k__Kingdom,p__Phylum,c__Class,o__Order,f__Family,g__Genus,s__Species,vitap_pi,vitap_confidence" > "\$OUT"

  python3 ${projectDir}/bin/fixvitapv0.4.py --vitap ${vitap_lineages} --header_map ${header_map} --out ${prefix}.vitap_evidence.csv --fallback none
  """
}

process RUN_GIANTHUNTER {

  tag "${prefix}"

  conda "${projectDir}/bin/gianthunter.yml"

  publishDir { "${params.outdir}/${prefix}/GiantHunter" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path(gianthunter_db)

  output:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path("${prefix}_gianthunter_prediction.tsv"),

  script:
  """
  set -euo pipefail

  outdir="${prefix}_gianthunter"
  mkdir -p "\$outdir"

  GiantHunter --contigs ${fasta} --dbdir ${gianthunter_db} --out "\$outdir"
  raw="./${prefix}_gianthunter/gianthunter"

  raw="\$outdir/final_prediction/gianthunter_prediction.tsv"

  # Always produce the file (even if empty) so downstream doesn't break
  if [[ -s "\$raw" ]]; then
    cp "\$raw" ${prefix}_gianthunter_prediction.tsv
  else
    # create empty-but-valid tsv with header
    echo -e "Accession\\tLength\\tGiantVirus\\tPotentialLineage\\tScore" > ${prefix}_gianthunter_prediction.tsv
  fi
 """
}

process PROCESS_GIANTHUNTER {

  tag "$prefix"

  conda 'bioconda::python=3.11 pandas'

  publishDir { "${params.outdir}/${prefix}/GiantHunter" }, mode: 'copy'

  input:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path(gianthunter_pred),


  output:
    tuple val(prefix),
          val("gianthunter"),
          path("${prefix}.gianthunter_evidence.csv")

  script:
  """
  set -euo pipefail

  OUT="${prefix}.gianthunter_evidence.csv"

  # Always create output with header so pipeline never breaks
  echo "seqid,d__Domain,r__Realm,k__Kingdom,p__Phylum,c__Class,o__Order,f__Family,g__Genus,s__Species,vitap_pi,vitap_confidence" > "\$OUT"

  awk -F "\t" '
  BEGIN {
    ranks[1]="d__"; ranks[2]="k__"; ranks[3]="p__";
    ranks[4]="c__"; ranks[5]="o__"; ranks[6]="f__";
    ranks[7]="g__"; ranks[8]="s__"
  }
  $4 ~ /Viruses/ {
    gsub(/superkingdom:/, "d__", $4)
    gsub(/kingdom:/, "k__", $4)
    gsub(/phylum:/, "p__", $4)
    gsub(/class:/, "c__", $4)
    gsub(/order:/, "o__", $4)
    gsub(/family:/, "f__", $4)
    gsub(/genus:/, "g__", $4)
    gsub(/species:/, "s__", $4)

    n = split($4, raw, ";")
    count = 0
    for (i = 1; i <= n; i++) {
        if (raw[i] ~ /^(d__|k__|p__|c__|o__|f__|g__|s__)/) {
            parts[++count] = raw[i]
        }
    }

    for (i = count+1; i <= 8; i++) {
        parts[i] = ranks[i] "unclassified"
    }

    lineage = parts[1]
    for (i = 2; i <= 8; i++) lineage = lineage ";" parts[i]

    print $1 "," lineage "," $5
  }
  '   ${gianthunter_pred} > "\$OUT"
  """
}

process RUN_VIRBOT {

  tag "$prefix"

  conda "${projectDir}/bin/virbot.yml"
  
  publishDir { "${params.outdir}/${prefix}/virbot" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path virbot_db

  output:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path("${prefix}.virbot.review"),
          path("${prefix}.pos_contig_score.csv")

  script:
  """
  set -euo pipefail

  outdir="${prefix}.virbot.review"
  mkdir -p "\$outdir"

  # Run virbot

  python VirBot.py -i ${norm_fasta} -o "\$outdir" -d ${virbot_db} [--sen] --threads ${params.threads}
 

  # virbot writes: <outdir>/pos_contig_score.csv
  raw="\$outdir/pos_contig_score.csv"

  # Always produce the file (even if empty) so downstream doesn't break
  if [[ -s "\$raw" ]]; then
    cp "\$raw" ${prefix}.pos_contig_score.csv
  else
    # create empty-but-valid tsv with header
    echo -e "Genome_ID\\tlineage\\tlineage_score/participation_index\\tConfidence_level" > ${prefix}.pos_contig_score.csv
  fi
  """
}

process RUN_CAT {

  tag "$prefix"

  conda "${projectDir}/bin/cat.yml"

  publishDir { "${params.outdir}/${prefix}/CAT" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path cat_db
    path cat_repo

  output:
    tuple val(prefix),
          val(type),
          path(norm_fasta),
          path(header_map),
          path(proteins_faa),
          path("${prefix}.cat_review"),
          path("${prefix}.cat_lineage.tsv")

  script:
  """
  set -euo pipefail

  outdir="${prefix}.cat_review"
  mkdir -p "\$outdir"

  CAT_pack contigs -i ${proteins_faa} -d ${cat_db} -o "/$outdir" -t ${params.threads}
  raw="\$outdir/contig_annotations.tsv"

  if [[ -s "\$raw" ]]; then
    cp "\$raw" ${prefix}.cat_lineage.tsv
  else
    # create empty-but-valid tsv with header
    echo -e "contig_id\\tlineage\\tscore\\tannotation_method\\tannotation_source\\tannotation_accession" > ${prefix}.cat_lineage.tsv
  fi
  """
}


workflow {

  /*
    Sample prefixes become part of sequence IDs and output filenames.
    Reject unsafe prefixes rather than silently changing them.
  */
  def validatePrefix = { rawPrefix, source ->
    if( rawPrefix == null || rawPrefix.toString().isEmpty() )
      error "Missing sample prefix (${source})."

    def prefix = rawPrefix.toString()
    def allowedPattern = /^[A-Za-z0-9][A-Za-z0-9._-]*$/

    if( !(prefix ==~ allowedPattern) ) {
      def invalidChars = prefix.replaceAll(/[A-Za-z0-9._-]/, '').toList().unique().join(' ')
      def suggestion = prefix.replaceAll(/[^A-Za-z0-9._-]/, '_')

      if( !(suggestion ==~ /^[A-Za-z0-9].*/) )
        suggestion = "sample${suggestion}"

      def detail = invalidChars
        ? " Invalid character(s): '${invalidChars}'."
        : " The first character must be a letter or number."

      error "Invalid sample prefix '${prefix}' (${source}).${detail} " +
            "Use only letters, numbers, periods, underscores, and hyphens; " +
            "the first character must be a letter or number. Suggested prefix: '${suggestion}'."
    }

    return prefix
  }

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

    def prefix = validatePrefix(params.prefix, '--prefix')
    ch_samples = Channel.of( tuple(prefix, t, fasta) )
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

        def prefix = row.prefix?.toString()
        def type   = row.type?.toString()?.trim()?.toLowerCase()
        def rel    = row.fasta?.toString()?.trim()

        if( !prefix || !type || !rel )
          error "CSV must have columns prefix,type,fasta with non-empty values. Bad row: ${row}"

        prefix = validatePrefix(prefix, "prefix_many row")

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

  // CHECK if ICTV list is present
  if( !file(params.ictv_csv).exists() )
      error "ICTV CSV not found: ${params.ictv_csv}"

  // CHECKV DB
  def CHECKV_DB_PATH = file("${params.dbdir}/checkv/checkv_db")

  Channel ch_checkv_db
  
  if( CHECKV_DB_PATH.exists() ) {
    ch_checkv_db = Channel.value(CHECKV_DB_PATH)
  } else {
    ch_checkv_db = CHECKV_DB().map { it -> CHECKV_DB_PATH }
  }

  // DEEP6 DB / SETUP
  def DEEP6_DB_PATH     = file(params.deep6_dir)
  def DEEP6_DB_SENTINEL = file("${params.deep6_dir}/.db_complete")

  Channel ch_deep6_db

  if( DEEP6_DB_SENTINEL.exists() ) {
    ch_deep6_db = Channel.value(DEEP6_DB_PATH)
  } else {
    ch_deep6_db = DEEP6_DB().map { it -> DEEP6_DB_PATH }
  }

  // GENOMAD DB
  def GENOMAD_DB_PATH     = file("${params.dbdir}/genomad/genomad_db")
  def GENOMAD_DB_SENTINEL = file("${params.dbdir}/genomad/genomad_db/.db_complete")

  Channel ch_genomad_db

  if( GENOMAD_DB_SENTINEL.exists() ) {
    ch_genomad_db = Channel.value(GENOMAD_DB_PATH)
  } else {
    ch_genomad_db = GENOMAD_DB().map { it -> GENOMAD_DB_PATH }
  }

  // VIRSORTER2 DB
  def VS2_DB_PATH     = file(params.vs2_dir)
  def VS2_DB_SENTINEL = file("${params.vs2_dir}/.db_complete")

  Channel ch_vs2_db

  if( VS2_DB_SENTINEL.exists() ) {
    ch_vs2_db = Channel.value(VS2_DB_PATH)
  } else {
    ch_vs2_db = VIRSORTER2_DB().map { it -> VS2_DB_PATH }
  }

  // CENOTE-TAKER3 DB
  def CT3_DB_PATH     = file(params.ct3_dir)
  def CT3_DB_SENTINEL = file("${params.ct3_dir}/.db_complete")

  Channel ch_ct3_db

  if( CT3_DB_SENTINEL.exists() ) {
    ch_ct3_db = Channel.value(CT3_DB_PATH)
  } else {
    ch_ct3_db = CENOTETAKER3_DB().map { it -> CT3_DB_PATH }
  }

  // VITAP DB 
  def VITAP_DB_PATH     = file("${params.vitap_dir}/DB_${params.vitap_db_label}")
  def VITAP_DB_SENTINEL = file("${params.vitap_dir}/DB_${params.vitap_db_label}/.db_complete")

  Channel ch_vitap_db

  if( VITAP_DB_SENTINEL.exists() ) {
    ch_vitap_db = Channel.value(VITAP_DB_PATH)
  } else {
    ch_vitap_db = VITAP_DB().map { it -> VITAP_DB_PATH }
  }

  // GIANTHUNTER DB
  def GIANTHUNTER_DB_PATH     = file("${params.gianthunter_dir}")
  def GIANTHUNTER_DB_SENTINEL = file("${params.gianthunter_dir}/.db_complete")

  Channel ch_gianthunter_db

  if( GIANTHUNTER_DB_SENTINEL.exists() ) {
    ch_gianthunter_db = Channel.value(GIANTHUNTER_DB_PATH)
  } else {
    ch_gianthunter_db = GIANTHUNTER_DB().map { it -> GIANTHUNTER_DB_PATH }
  } 

  // VirBot  repo and database path check
  def VIRBOT_REPO_PATH = file(params.virbot_repo_dir)
  if( !VIRBOT_REPO_PATH.exists() ) {
    error "VirBot repo directory not found: ${params.virbot_repo_dir}"
  }

  def VIRBOT_DB_PATH = file(params.virbot_db_dir)
  if( !VIRBOT_DB_PATH.exists() ) {
    error "VirBot database directory not found: ${params.virbot_db_dir}"
  }

  Channel ch_virbot_repo = Channel.value(VIRBOT_REPO_PATH)
  Channel ch_virbot_db   = Channel.value(VIRBOT_DB_PATH)

  // VCONTACT3 DB
  def VCONTACT3_DB_PATH     = file(params.vcontact3_db_dir)
  def VCONTACT3_DB_SENTINEL = file("${params.vcontact3_db_dir}/.db_complete")

  Channel ch_vcontact3_db

  if( VCONTACT3_DB_SENTINEL.exists() ) {
    ch_vcontact3_db = Channel.value(VCONTACT3_DB_PATH)
  } else {
    ch_vcontact3_db = VCONTACT3_DB().map { it -> VCONTACT3_DB_PATH }
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
  ch_deep6_ev = RUN_DEEP6(branched.rna, ch_deep6_db)

  // GENOMAD PROCESSES
  ch_genomad_raw = RUN_GENOMAD(ch_orf, ch_genomad_db)
  ch_genomad_raw.view { p,t,nf,hm,faa,review,vs,ps ->
    "GENOMAD sample=${p}\tvirus_summary=${vs.name}\tplasmid_summary=${ps.name}"
  }
  ch_genomad_ev = GENOMAD_PROCESSING(ch_genomad_raw, ictv_csv_ch)

  // VIRSORTER2 PROCESSES
  ch_vs2_ev = RUN_VIRSORTER2(ch_orf, ch_vs2_db)

  // CENOTETAKER3 PROCESSES
  ch_ct3_raw = RUN_CENOTETAKER3(ch_orf, ch_ct3_db)
  ch_ct3_ev  = PROCESS_CENOTETAKER3(ch_ct3_raw)

  // VITAP PROCESSES
  ch_vitap_raw = RUN_VITAP(ch_orf, ch_vitap_db)
  ch_vitap_ev  = PROCESS_VITAP(ch_vitap_raw, ictv_csv_ch)

  // GIANTHUNTER PROCESSES
  ch_gianthunter_raw = RUN_GIANTHUNTER(ch_orf, ch_gianthunter_db)
  ch_gianthunter_ev = PROCESS_GIANTHUNTER(ch_gianthunter_raw)

  // VIRBOT PROCESSSES


// BELOW IS THE TEMPLATE TO MERGE EVIDENCE FROM ALL PROGRAMS RAN ON FASTA FILE

// core tuple from ORF stage: (prefix, type, norm_fasta, header_map, proteins_faa)
ch_core = ch_orf.map { p, t, nf, hm, faa -> tuple(p, t, nf, hm, faa) }

// gather all evidence records (prefix, tool, file)
ch_all_evidence = ch_genomad_ev
  .mix(ch_checkv_ev)
  .mix(ch_deep6_ev)
  .mix(ch_vs2_ev)
  .mix(ch_ct3_ev)
  .mix(ch_vitap_ev)
  .mix(ch_gianthunter_ev)
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
