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
    tuple val(prefix),
          val(type),
          path("${prefix}.normalized.fasta"),
          path("${prefix}.header_map.tsv")

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

process RUN_GENOMAD {

  tag "$prefix"

  conda 'bioconda::genomad'

  publishDir { "${params.outdir}/${prefix}/genomad" }, mode: 'copy'

  input:
    tuple val(prefix), val(type), path(norm_fasta), path(header_map), path(proteins_faa)
    path genomad_db

  output:
    tuple val(prefix), path("${prefix}.genomad.evidence.tsv")
    path("${prefix}.genomad.review")

  script:
  """
  mkdir -p ${prefix}.genomad.review

  # Example invocation (you’ll fill in your real options)
  # genomad end-to-end --cleanup --threads ${params.threads} ${norm_fasta} ${prefix}.genomad.review ${genomad_db}

  # Placeholder evidence file for now
  echo -e "contig\\ttool\\tscore" > ${prefix}.genomad.evidence.tsv
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

// stable final location you want to use in all runs
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



}
