#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=gff_compare
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL


#1.lyric vs isoq
#2. ref vs lyric
#3. ref vs isoq

export PATH=/users/rg/baygun/gff_compare/gffcompare:$PATH

gffcompare -r lyric.gtf -o lyric_vs_isoquant isoquant_transcript_models.gtf

gffcompare -r *.annotation.gtf -o ref_lyric_gff_compare lyric.gtf

gffcompare -r *.annotation.gtf -o ref_isoquant_gff_compare isoquant_transcript_models.gtf

