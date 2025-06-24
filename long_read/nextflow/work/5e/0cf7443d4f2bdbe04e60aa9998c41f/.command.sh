#!/bin/bash -ue
mkdir -p "result_nanoplot/FAX00328_20240115_1815_946aeace"
NanoPlot --ubam "FAX00328_20240115_1815_946aeace.bam" -o result_nanoplot/FAX00328_20240115_1815_946aeace --N50 -threads 4
