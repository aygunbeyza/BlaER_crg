#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=tss
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

module load BEDTools/2.31.1

# 1. TSS 
zcat 120_lyric.gtf.gz | \
awk 'BEGIN{OFS="\t"} {
  if ($7 == "+" && $4 > 52) print $1, $4 - 1, $4, ".", ".", "+", "t120";
  else if ($7 == "-" && $5 > 52) print $1, $5 - 1, $5, ".", ".", "-", "t120";
}' > tss_120.bed

# 1. TSS 
zcat 0_lyric.gtf.gz | \
awk 'BEGIN{OFS="\t"} {
  if ($7 == "+" && $4 > 52) print $1, $4 - 1, $4, ".", ".", "+", "t0";
  else if ($7 == "-" && $5 > 52) print $1, $5 - 1, $5, ".", ".", "-", "t0";
}' > tss_0.bed

# 2. ±50bp  
bedtools slop -i /tss_120.bed -g /hg38.genome -b 50 > /tss120_50bp.bed

# 2. ±50bp
bedtools slop -i /tss_0.bed -g /hg38.genome -b 50 > /tss0_50bp.bed

# 3. CAGE peaks ile kesişim (t120 için)
bedtools intersect -a /tss120_50bp.bed \
-b CAGE_peaks_phase1and2.bed \
-u > supported_tss_t120.bed

# 3. CAGE peaks ile kesişim (t0 için)
bedtools intersect -a /tss0_50bp.bed \
-b /CAGE_peaks_phase1and2.bed \
-u > supported_tss_t0.bed

# 1. Desteklenmeyen TSS'ler - t120 için (CAGE peaks ile örtüşmeyenler)
bedtools intersect -a /tss120_50bp.bed \
-b /CAGE_peaks_phase1and2.bed \
-v > unsupported_tss_t120.bed

# 2. Desteklenmeyen TSS'ler - t0 için (CAGE peaks ile örtüşmeyenler)
bedtools intersect -a /tss0_50bp.bed \
-b /CAGE_peaks_phase1and2.bed \
-v > unsupported_tss_t0.bed

module load matplotlib/3.8.2-gfbf-2023b
python3 tss_graph.py
