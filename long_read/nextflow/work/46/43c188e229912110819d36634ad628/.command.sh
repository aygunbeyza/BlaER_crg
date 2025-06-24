#!/bin/bash -ue
mkdir -p "result_nanoplot/FAX02360_20240115_1630_e5122b92"
NanoPlot --ubam "FAX02360_20240115_1630_e5122b92.bam" -o result_nanoplot/FAX02360_20240115_1630_e5122b92 --N50 --threads 4
