#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=96:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=isoquant
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

conda init
source ~/.bashrc
conda activate isoquant

isoquant.py --reference /users/rg/baygun/BlaER_crg/ref/anotation_and_hg/hg38.sorted.fa.gz \
--genedb /users/rg/baygun/BlaER_crg/ref/anotation_and_hg/gencode.v47.primary_assembly.annotation.gtf.gz \
--complete_genedb \
--bam /users/rg/baygun/BlaER_crg/files/long_read_files/align_bam_files_lyric/ont-Crg-sc_HpreCap_0+_BLaER101Rep1.bam \
--data_type nanopore -o /users/rg/baygun/BlaER_crg/long_read/IsoQuant/result_isoquant_t0 \
--thread 8


conda deactivate

