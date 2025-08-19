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

module load STAR/2.7.11b-GCC-12.3.0


STAR --runThreadN 6 \
        --genomeDir "/users/rg/baygun/BlaER_crg/short_read/star/deneme_genome_index" \
        --readFilesIn /users/rg/baygun/BlaER_crg/files/short_read_files/scH000_lib_02535AAF_ATGACGTCGC-ATCCTGACCT_R2_001.fastq /users/rg/baygun/BlaER_crg/files/short_read_files/scH000_lib_02535AAF_ATGACGTCGC-ATCCTGACCT_R1_001.fastq \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix deneme_t0 \
        --soloType Droplet \
        --soloCBwhitelist whitelist.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 12 \
	--soloFeatures Gene GeneFull \
        --outSAMattributes All
