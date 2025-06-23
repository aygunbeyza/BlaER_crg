#!/bin/bash -ue
mkdir -p $file_id
NanoPlot --ubam FAX02360_20240115_1630_e5122b92.bam -o $file_id --threads 4 --NS0
mv $file_id/* .
