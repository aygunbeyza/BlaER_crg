#!/bin/bash -ue
mkdir -p $file_id
NanoPlot --ubam FAX00264_20240115_1817_5ed29589.bam -o $file_id --threads 4 --NS0
mv $file_id/* .
