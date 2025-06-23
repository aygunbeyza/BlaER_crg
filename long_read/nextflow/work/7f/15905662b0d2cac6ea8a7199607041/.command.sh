#!/bin/bash -ue
file_id=`basename FAX00264_20240115_1817_5ed29589.bam .bam`
mkdir -p result_nanoplot/${file_id}
NanoPlot --ubam FAX00264_20240115_1817_5ed29589.bam -o result_nanoplot/${file_id} --threads 4 --N50
