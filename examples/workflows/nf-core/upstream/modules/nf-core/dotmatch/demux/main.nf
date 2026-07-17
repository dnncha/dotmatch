process DOTMATCH_DEMUX {
    tag "$meta.id"
    label 'process_medium'

    cpus   { task.ext.cpus   ?: 4 }
    memory { task.ext.memory ?: 8.GB }
    time   { task.ext.time   ?: 4.h  }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dotmatch:0.1.9--h*' :
        'biocontainers/dotmatch:0.1.9--h*' }"

    input:
    tuple val(meta), path(reads), path(barcodes)
    val barcode_start
    val barcode_length
    val k
    val metric

    output:
    tuple val(meta), path("demuxed"), emit: demuxed
    tuple val(meta), path("*.summary.json"), emit: summary
    tuple val(meta), path("*.assignments.tsv"), emit: assignments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    dotmatch demux \\
      --barcodes ${barcodes} \\
      --reads ${reads} \\
      --barcode-start ${barcode_start} \\
      --barcode-length ${barcode_length} \\
      --k ${k} \\
      --metric ${metric} \\
      --ambiguity-policy radius \\
      --out-dir demuxed \\
      --summary ${prefix}.summary.json \\
      --assignments ${prefix}.assignments.tsv \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
