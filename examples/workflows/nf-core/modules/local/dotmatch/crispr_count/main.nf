process DOTMATCH_CRISPR_COUNT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads), path(library)
    val guide_start
    val guide_length
    val k
    val metric

    output:
    tuple val(meta), path("*.counts.mageck.tsv"), emit: counts
    tuple val(meta), path("*.summary.json"), emit: summary
    tuple val(meta), path("*.sample_qc.tsv"), emit: sample_qc
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    printf 'sample_id\\tfastq\\n' > samples.tsv
    printf '%s\\t%s\\n' '${meta.id}' '${reads}' >> samples.tsv

    dotmatch crispr-count \\
      --library ${library} \\
      --samples samples.tsv \\
      --guide-start ${guide_start} \\
      --guide-length ${guide_length} \\
      --k ${k} \\
      --metric ${metric} \\
      --out ${prefix}.counts.mageck.tsv \\
      --summary ${prefix}.summary.json \\
      --sample-qc ${prefix}.sample_qc.tsv \\
      --ambiguous discard \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
