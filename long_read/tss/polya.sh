#!/bin/bash
#SBATCH --output=polya.%A.%a.out
#SBATCH --error=polya.%A.%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=polya
#SBATCH --mail-user=reyza.abeygun@crg.eu
#SBATCH --mail-type=FAIL

module load BEDTools/2.31.1

zcat t120_lyric.gtf.gz | \
awk 'BEGIN{OFS="\t"} { if ($7 == "+" && $5 > 52) print $1, $5 - 1, $5, $10, ".", $7, "p120"; else if ($7 == "-" && $4 > 52) print $1, $4 - 1, $4, $10, ".", $7, "p120"; }' > polya_p120.bed

zcat t0_lyric.gtfgz | \
awk 'BEGIN{OFS="\t"} { if ($7 == "+" && $5 > 52) print $1, $5 - 1, $5, $10, ".", $7, "p0"; else if ($7 == "-" && $4 > 52) print $1, $4 - 1, $4, $10, ".", $7, "p0"; }' > polya_p0.bed

#50bp (p120)
bedtools slop -i polya_p120.bed -g /hg38.genome -b 50 > polya_p120_50bp.bed

bedtools slop -i polya_p0.bed -g /hg38.genome -b 50 > polya_p0_50bp.bed

#CAGE peaks(p120)
bedtools intersect -a polya_p120_50bp.bed \
  -b /CAGE_peaks_phase1and2.bed \
  > supported_polya_p120.bed

bedtools intersect -a polya_p0_50bp.bed \
  -b /CAGE_peaks_phase1and2.bed \
  > supported_polya_p0.bed

bedtools intersect -v -a polya_p120_50bp.bed \
  -b /CAGE_peaks_phase1and2.bed \
  > unsupported_polya_p120.bed

bedtools intersect -v -a polya_p0_50bp.bed \
  -b /CAGE_peaks_phase1and2.bed \
  > unsupported_polya_p0.bed

module load matplotlib/3.8.2-gfbf-2023b
python3 polya_graph.py
