#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// parameters
params.bam_root       = "/users/rg/baygun/BlaER_crg/files/long_read_files"
params.out_nanoplot   = "/users/rg/baygun/BlaER_crg/long_read/nextflow/result_nanoplot"
params.out_epi2me     = "/users/rg/baygun/BlaER_crg/long_read/nextflow/results_epi2me"
params.ref_genome_dir = "/users/rg/baygun/BlaER_crg/ref"

process NanoPlot {
    input:
    path input_bam

    script:
    """
    mkdir -p "${params.out_nanoplot}/${input_bam.simpleName}"
    NanoPlot --ubam "${input_bam}" -o "${params.out_nanoplot}/${input_bam.simpleName}" --N50 --threads 4
    """
}

process Epi2me {
    tag "${sample_dir.simpleName}"

    input:
    path sample_dir

    script:
    """
    nextflow run epi2me-labs/wf-single-cell -r v3.3.0 \\
        --expected_cells 10000 \\
        --bam ${sample_dir} \\
        --kit '3prime:v3' \\
        --ref_genome_dir ${params.ref_genome_dir} \\
        --out_dir ${params.out_epi2me}/${sample_dir.simpleName} \\
        -profile standard,singularity \\
        -without-docker
    """
}

workflow {
    // NanoPlot: tüm bam dosyaları (t0 & t120)
    Channel
        .fromPath("${params.bam_root}/**/*.bam")
        | NanoPlot

    // epi2me: her raw klasörü bir sample gibi
    Channel
        .fromPath("${params.bam_root}/*_raw_files")
        .filter { it.isDirectory() }
        | Epi2me
}

