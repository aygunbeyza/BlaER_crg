#!/bin/bash -ue
file_id=`basename FAX02360_20240115_1630_e5122b92.bam .bam`
mkdir -p result_nanoplot/${file_id}
NanoPlot --ubam FAX02360_20240115_1630_e5122b92.bam -o result_nanoplot/${file_id} --threads 4 --N50
