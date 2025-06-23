#!/bin/bash -ue
mkdir -p FAX00328_20240115_1815_946aeace
NanoPlot --ubam FAX00328_20240115_1815_946aeace.bam -o FAX00328_20240115_1815_946aeace --threads 4 --NS0
mv FAX00328_20240115_1815_946aeace/* .
