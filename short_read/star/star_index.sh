#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=deneme_star
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

module load STAR/2.7.11b-GCC-12.3.0

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir "/genome_index/" \
--genomeFastaFiles "/hg38.fa" \
--sjdbGTFfile "/.annotation.gtf" \
--sjdbOverhang 89
