workflow {
    Channel.fromPath("${params.bam_folder}/*.bam")
        .set { bam_files }
    run_nanoplot(bam_files)
}

process run_nanoplot {
    publishDir "result_nanoplot", mode: 'copy'
    tag "$bam_file"

    input:
    path bam_file

    output:
    path "*"

    script:
    def file_id = bam_file.getBaseName()  // Properly defining file_id
    """
    mkdir -p $file_id
    NanoPlot --ubam $bam_file -o $file_id --threads 4 --NS0
    mv $file_id/* .
    """
}

