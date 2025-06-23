nextflow.enable.dsl=2

workflow {
  Channel.fromPath("${params.bam_folder}/*.bam")
         .set { bam_files }

  run_nanoplot(bam_files)
}

process run_nanoplot {
  tag "$bam_file"

  input:
    path bam_file

  output:
    path "result_nanoplot/${bam_file.simpleName}/*", optional: true

  script:
    """
    file_id=`basename $bam_file .bam`
    mkdir -p result_nanoplot/\${file_id}
    NanoPlot --ubam $bam_file -o result_nanoplot/\${file_id} --threads 4 --N50
    """
}

