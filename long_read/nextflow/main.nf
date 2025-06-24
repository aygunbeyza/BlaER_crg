#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process NanoPlot {
    input:
    path input_bam
    
    script:
    """
    mkdir -p "${input_bam.baseName}"
    NanoPlot --ubam "${input_bam}" --N50 --threads 4
    """
}


workflow {
    Channel.fromPath("/users/rg/baygun/BlaER_crg/files/long_read_files/*.bam") | NanoPlot
}
