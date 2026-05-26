process DOTMATCH_ASSAY_RUN {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assay_spec), path(assay_inputs)

    output:
    tuple val(meta), path("assay_report.html"), emit: assay_report
    tuple val(meta), path("assay_manifest.json"), emit: assay_manifest
    tuple val(meta), path("assay_manifest.summary.tsv"), emit: assay_manifest_summary
    tuple val(meta), path("sample_qc.tsv"), emit: sample_qc
    tuple val(meta), path("crispr_qc.html"), emit: crispr_qc_report
    tuple val(meta), path("crispr_qc.json"), emit: crispr_qc_json
    tuple val(meta), path("crispr_qc.summary.tsv"), emit: crispr_qc_summary
    tuple val(meta), path("counts.mageck.tsv"), emit: counts
    tuple val(meta), path("summary.json"), emit: summary
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    for input_file in ${assay_inputs}; do
      ln -sf "\${input_file}" "\$(basename "\${input_file}")"
    done

    dotmatch assay run ${assay_spec} ${args}

    cp assay_out/assay_report.html assay_report.html
    cp assay_out/assay_manifest.json assay_manifest.json
    cp assay_out/assay_manifest.summary.tsv assay_manifest.summary.tsv
    cp assay_out/sample_qc.tsv sample_qc.tsv
    cp assay_out/crispr_qc.html crispr_qc.html
    cp assay_out/crispr_qc.json crispr_qc.json
    cp assay_out/crispr_qc.summary.tsv crispr_qc.summary.tsv
    cp assay_out/counts.mageck.tsv counts.mageck.tsv
    cp assay_out/summary.json summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      dotmatch: "\$(dotmatch --version | sed 's/^dotmatch //')"
    END_VERSIONS
    """
}
