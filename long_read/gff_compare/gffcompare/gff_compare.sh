#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=gff_compare
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL
export PATH=/users/rg/baygun/gff_compare/gffcompare:$PATH

zcat /users/rg/baygun/BlaER_crg/files/long_read_files/lyric_gtf/ont-Crg-sc_Hv3_0+_BLaER101Rep1.HiSS.tmerge.min2reads.splicing_status-all.endSupport-all.gff.gz > lyric.gtf

#gffcompare -r /users/rg/baygun/BlaER_crg/long_read/sqanti/ref/gencode.v47.primary_assembly.annotation.gtf -o t120_gff_compare lyric.gtf /users/rg/baygun/BlaER_crg/long_read/IsoQuant/result_isoquant_t120/OUT/OUT.transcript_models.gtf

gffcompare -r lyric.gtf -o lyric_vs_isoquant isoquant_transcript_models.gtf

#gffcompare -r /users/rg/baygun/BlaER_crg/long_read/sqanti/ref/gencode.v47.primary_assembly.annotation.gtf -o ref_lyric_gff_compare /users/rg/baygun/BlaER_crg/long_read/gff_compare/gffcompare/lyric.gtf

#gffcompare -r /users/rg/baygun/BlaER_crg/long_read/sqanti/ref/gencode.v47.primary_assembly.annotation.gtf -o ref_isoquant_gff_compare isoquant_transcript_models.gtf
