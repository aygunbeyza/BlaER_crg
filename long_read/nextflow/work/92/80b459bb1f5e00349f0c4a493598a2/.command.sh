#!/bin/bash -ue
file_id=`basename FAX02366_20240115_1632_794829ec.bam .bam`
mkdir -p result_nanoplot/${file_id}
NanoPlot --ubam FAX02366_20240115_1632_794829ec.bam -o result_nanoplot/${file_id} --threads 4 --N50
