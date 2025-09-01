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

isoquant.py --reference /hg38.sorted.fa.gz \
--genedb ssembly.annotation.gtf.gz \
--complete_genedb \
--bam lyric.bam \
--data_type nanopore -o /result_isoquant_t0 \
--thread 8


conda deactivate

