nextflow.enable.dsl=2

sample_rows = []
file(params.samples).readLines().eachWithIndex { line, row_index ->
    def trimmed = line.trim()
    if (trimmed && !trimmed.startsWith('#') && row_index > 0) {
        def fields = line.split('\t')
        if (fields.size() < 2) {
            throw new IllegalArgumentException("samples.tsv must contain sample_id and fastq columns")
        }
        sample_rows << [sample_id: fields[0], fastq: file(fields[1])]
    }
}

if (sample_rows.isEmpty()) {
    throw new IllegalArgumentException("samples.tsv must contain at least one sample")
}

process DOTMATCH_CRISPR_COUNT {
    tag "dotmatch_crispr_count"
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path library
    val sample_ids
    path fastqs

    output:
    path "counts.mageck.tsv", emit: counts
    path "summary.json", emit: summary
    path "sample_qc.tsv", emit: sample_qc

    script:
    def indelWindow = params.indel_window as int
    def indelArg = params.metric == 'levenshtein' && indelWindow > 0 ? "--indel-window ${indelWindow}" : ""
    def sampleSheetRows = sample_ids.indices.collect { i -> "${sample_ids[i]}\t${fastqs[i].name}" }.join('\n')
    """
    printf 'sample_id\\tfastq\\n' > samples.tsv
    printf '%s\\n' '${sampleSheetRows}' >> samples.tsv

    dotmatch crispr-count \
      --library ${library} \
      --samples samples.tsv \
      --guide-start ${params.guide_start} \
      --guide-length ${params.guide_length} \
      --k ${params.k} \
      --metric ${params.metric} \
      ${indelArg} \
      --out counts.mageck.tsv \
      --summary summary.json \
      --sample-qc sample_qc.tsv \
      --ambiguous discard
    """
}

process DOTMATCH_ASSAY_RUN {
    tag "dotmatch_assay_run"
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path library
    val sample_ids
    path fastqs

    output:
    path "assay_report.html", emit: assay_report
    path "assay_manifest.json", emit: assay_manifest
    path "assay_manifest.summary.tsv", emit: assay_manifest_summary
    path "sample_qc.tsv", emit: assay_sample_qc
    path "counts.mageck.tsv", emit: assay_counts
    path "summary.json", emit: assay_summary

    script:
    def indelWindow = params.indel_window as int
    def indelLine = params.metric == 'levenshtein' && indelWindow > 0 ? "indel_window = ${indelWindow}" : ""
    def sampleBlocks = sample_ids.indices.collect { i -> """
[[samples]]
id = "${sample_ids[i]}"
fastq = "${fastqs[i].name}"
""" }.join('\n')
    """
    cat > assay.toml <<'ASSAY'
    schema_version = 1
    status = "ready"
    mode = "count"
    assay_type = "crispr"
    targets = "${library}"

    ${sampleBlocks}
    [run]
    out_dir = "assay_out"
    threads = ${task.cpus}

    [extract]
    start = ${params.guide_start}
    length = ${params.guide_length}

    [assignment]
    k = ${params.k}
    metric = "${params.metric}"
    ambiguous = "discard"
    ${indelLine}

    [outputs]
    format = "mageck"
    ASSAY

    dotmatch assay run assay.toml
    cp assay_out/assay_report.html assay_report.html
    cp assay_out/assay_manifest.json assay_manifest.json
    cp assay_out/assay_manifest.summary.tsv assay_manifest.summary.tsv
    cp assay_out/sample_qc.tsv sample_qc.tsv
    cp assay_out/counts.mageck.tsv counts.mageck.tsv
    cp assay_out/summary.json summary.json
    """
}

workflow {
    sample_ids_ch = Channel.value(sample_rows.collect { it.sample_id })
    fastqs_ch = Channel.value(sample_rows.collect { it.fastq })
    library_ch = Channel.fromPath(params.library, checkIfExists: true)

    DOTMATCH_CRISPR_COUNT(library_ch, sample_ids_ch, fastqs_ch)
    DOTMATCH_ASSAY_RUN(library_ch, sample_ids_ch, fastqs_ch)
}
