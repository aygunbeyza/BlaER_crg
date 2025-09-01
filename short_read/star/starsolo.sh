#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --mem=47G
#SBATCH --cpus-per-task=6
#SBATCH --job-name=star_solo
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

#        --readFilesCommand zcat \
#This script sets the genome directory, creates an index folder, and runs STAR to generate a genome index for the human genome (hg38) using the FASTA file and GTF annotation file.
#The command includes a commented line for setting splice junction overhangs, which can be used to improve alignment accuracy.
#NOTE FOR STARSOLO: for t0 barcode are umi in r2 be careful
#t0: R2-R1
#t120: R1 -R2 but when i use like that ı had a error that your barcode and umi file length is 90 so i wrote firstly r2

module load STAR/2.7.11b-GCC-12.3.0


STAR --runThreadN 6 \
        --genomeDir "/genome_index" \
        --readFilesIn /R2_001.fastq /_R1_001.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix t0_with_genefull_ \
        --soloType Droplet \
        --soloCBwhitelist whitelist.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloFeatures Gene GeneFull \
        --outSAMattributes All


STAR --runThreadN 6 \
        --genomeDir "/genome_index" \
        --readFilesIn /R2_001.fastq  /R1_001.fastq 
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix t120_with_genefull_ \
        --soloType Droplet \
        --soloCBwhitelist whitelist.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
        --soloFeatures Gene GeneFull \
        --outSAMattributes All
starsolo.sh

