#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=nanoplot_4
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

conda init
conda activate qc3
mkdir /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX00264_20240115_1817_5ed29589
mkdir /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX00328_20240115_1815_946aeace
mkdir /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX02360_20240115_1630_e5122b92
mkdir /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX02366_20240115_1632_794829ec

NanoPlot --ubam /users/rg/baygun/BlaER/files/long_read_files/FAX00264_20240115_1817_5ed29589.bam \
-o /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX00264_20240115_1817_5ed29589 \
--N50 \
--threads 4

NanoPlot --ubam /users/rg/baygun/BlaER/files/long_read_files/FAX00328_20240115_1815_946aeace.bam \
-o /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX00328_20240115_1815_946aeace \
--N50 \
--threads 4

NanoPlot --ubam /users/rg/baygun/BlaER/files/long_read_files/FAX02360_20240115_1630_e5122b92.bam \
-o /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX02360_20240115_1630_e5122b92 \
--N50 \
--threads 4

NanoPlot --ubam /users/rg/baygun/BlaER/files/long_read_files/FAX02366_20240115_1632_794829ec.bam \
-o /users/rg/baygun/BlaER/long_read/nanoplot/result_nanoplot/FAX02366_20240115_1632_794829ec \
--N50 \
--threads 4







conda deactivate
