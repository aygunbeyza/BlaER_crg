#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process NanoPlot {
    input:
    path input_bam

    output:
    path "result_nanoplot/${input_bam.baseName}"

    script:
    """
    mkdir -p "result_nanoplot/${input_bam.baseName}"
    NanoPlot --ubam "${input_bam}" -o result_nanoplot/${input_bam.baseName} --N50 --threads 4
    """
}


workflow {
    Channel.fromPath("/users/rg/baygun/BlaER_crg/files/long_read_files/*.bam") | NanoPlot
}
