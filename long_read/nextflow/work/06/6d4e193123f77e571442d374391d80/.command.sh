#!/bin/bash -ue
file_id=`basename FAX00328_20240115_1815_946aeace.bam .bam`
mkdir -p result_nanoplot/${file_id}
NanoPlot --ubam FAX00328_20240115_1815_946aeace.bam -o result_nanoplot/${file_id} --threads 4 --N50
