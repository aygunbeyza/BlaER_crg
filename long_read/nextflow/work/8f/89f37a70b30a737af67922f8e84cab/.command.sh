#!/bin/bash -ue
mkdir -p "result_nanoplot/FAX00264_20240115_1817_5ed29589"
NanoPlot --ubam "FAX00264_20240115_1817_5ed29589.bam" -o result_nanoplot/FAX00264_20240115_1817_5ed29589 --N50 --threads 4
