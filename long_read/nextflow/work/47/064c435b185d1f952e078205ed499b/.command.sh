#!/bin/bash -ue
mkdir -p $file_id
NanoPlot --ubam FAX02366_20240115_1632_794829ec.bam -o $file_id --threads 4 --NS0
mv $file_id/* .
