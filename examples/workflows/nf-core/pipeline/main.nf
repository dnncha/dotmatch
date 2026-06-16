#!/usr/bin/env nextflow

/*
 * Example end-to-end nf-core-style pipeline for DotMatch
 *
 * Demonstrates both the crispr_count and assay_run modules on the shared fixtures.
 * This is a minimal working example for testing/adoption.
 *
 * Usage:
 *   nextflow run . --outdir results
 *
 * In a real pipeline, use a samplesheet and proper input channels.
 */

nextflow.enable.dsl = 2

// Include the local module candidates (relative from pipeline/ dir)
include { DOTMATCH_CRISPR_COUNT } from '../modules/local/dotmatch/crispr_count/main'
include { DOTMATCH_ASSAY_RUN }    from '../modules/local/dotmatch/assay_run/main'

workflow {
    // Demo for CRISPR count module using fixture
    ch_crispr_input = Channel.of([
        [ id: 'demo_crispr' ],
        file('examples/workflows/fixtures/sample_a.fastq'),
        file('examples/workflows/fixtures/crispr_library.csv')
    ])

    DOTMATCH_CRISPR_COUNT (
        ch_crispr_input,
        0,   // guide_start
        4,   // guide_length (matches fixture)
        1,   // k
        'hamming'
    )

    // Demo for full AssaySpec run
    ch_assay_spec = Channel.of([
        [ id: 'demo_assay' ],
        file('examples/workflows/fixtures/crispr_assay.toml'),
        [
            file('examples/workflows/fixtures/crispr_library.csv'),
            file('examples/workflows/fixtures/sample_a.fastq'),
            file('examples/workflows/fixtures/sample_b.fastq')
        ]
    ])

    DOTMATCH_ASSAY_RUN ( ch_assay_spec )

    // Emit for inspection
    DOTMATCH_CRISPR_COUNT.out.counts.view { "CRISPR counts: $it" }
    DOTMATCH_ASSAY_RUN.out.assay_report.view { "Assay report: $it" }
}

workflow.onComplete {
    println "DotMatch example pipeline complete. Results in ${params.outdir ?: 'work'}"
}