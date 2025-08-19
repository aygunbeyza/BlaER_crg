#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=16:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=nanoplot_4
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

conda init
source ~/.bashrc
conda activate qc3

IN_DIR="/users/rg/baygun/BlaER_crg/files/long_read_files/all_raw_bams"
OUT_DIR="/users/rg/baygun/BlaER_crg/long_read/nanoplot/result_nanoplot"

mkdir -p "$OUT_DIR"

# for all bam files
for bam in "$IN_DIR"/*.bam; do
    base=$(basename "$bam" .bam)
    out="$OUT_DIR/$base"
    mkdir -p "$out"

    NanoPlot --ubam "$bam" \
             -o "$out" \
             --N50 \
             --threads 4
done

conda deactivate
