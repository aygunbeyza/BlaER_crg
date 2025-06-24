#!/bin/bash -ue
mkdir -p "result_nanoplot/FAX02366_20240115_1632_794829ec"
NanoPlot --ubam "FAX02366_20240115_1632_794829ec.bam" -o result_nanoplot/FAX02366_20240115_1632_794829ec --N50 --threads 4
