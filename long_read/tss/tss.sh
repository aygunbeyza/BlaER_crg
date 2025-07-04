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

# 1. TSS çıkarma - İlk veri seti (t120)
zcat /users/rg/baygun/BlaER_crg/files/long_read_files/lyric_gtf/ont-Crg-sc_Hv3_0+_BLaER101Rep1.HiSS.tmerge.min2reads.splicing_status-all.endSupport-all.gff.gz | \
awk 'BEGIN{OFS="\t"} {
  if ($7 == "+" && $4 > 52) print $1, $4 - 1, $4, ".", ".", "+", "t120";
  else if ($7 == "-" && $5 > 52) print $1, $5 - 1, $5, ".", ".", "-", "t120";
}' > tss_120.bed

# 1. TSS çıkarma - İkinci veri seti (t0)
zcat /users/rg/baygun/BlaER_crg/files/long_read_files/lyric_gtf/ont-Crg-sc_HpreCap_0+_BLaER101Rep1.HiSS.tmerge.min2reads.splicing_status-all.endSupport-all.gff.gz | \
awk 'BEGIN{OFS="\t"} {
  if ($7 == "+" && $4 > 52) print $1, $4 - 1, $4, ".", ".", "+", "t0";
  else if ($7 == "-" && $5 > 52) print $1, $5 - 1, $5, ".", ".", "-", "t0";
}' > tss_0.bed

# 2. ±50bp genişletme (t120 için)
bedtools slop -i /users/rg/baygun/BlaER_crg/long_read/tss/tss_120.bed -g /users/rg/baygun/BlaER_crg/long_read/tss/hg38.genome -b 50 > /users/rg/baygun/BlaER_crg/long_read/tss/tss120_50bp.bed

# 2. ±50bp genişletme (t0 için)
bedtools slop -i /users/rg/baygun/BlaER_crg/long_read/tss/tss_0.bed -g /users/rg/baygun/BlaER_crg/long_read/tss/hg38.genome -b 50 > /users/rg/baygun/BlaER_crg/long_read/tss/tss0_50bp.bed

# 3. CAGE peaks ile kesişim (t120 için)
bedtools intersect -a /users/rg/baygun/BlaER_crg/long_read/tss/tss120_50bp.bed \
-b /users/project/gencode_006070_no_backup/jlagarde/fantom5/hg38_fair+new_CAGE_peaks_phase1and2.bed \
-u > supported_tss_t120.bed

# 3. CAGE peaks ile kesişim (t0 için)
bedtools intersect -a /users/rg/baygun/BlaER_crg/long_read/tss/tss0_50bp.bed \
-b /users/project/gencode_006070_no_backup/jlagarde/fantom5/hg38_fair+new_CAGE_peaks_phase1and2.bed \
-u > supported_tss_t0.bed

# 1. Desteklenmeyen TSS'ler - t120 için (CAGE peaks ile örtüşmeyenler)
bedtools intersect -a /users/rg/baygun/BlaER_crg/long_read/tss/tss120_50bp.bed \
-b /users/project/gencode_006070_no_backup/jlagarde/fantom5/hg38_fair+new_CAGE_peaks_phase1and2.bed \
-v > unsupported_tss_t120.bed

# 2. Desteklenmeyen TSS'ler - t0 için (CAGE peaks ile örtüşmeyenler)
bedtools intersect -a /users/rg/baygun/BlaER_crg/long_read/tss/tss0_50bp.bed \
-b /users/project/gencode_006070_no_backup/jlagarde/fantom5/hg38_fair+new_CAGE_peaks_phase1and2.bed \
-v > unsupported_tss_t0.bed

module load matplotlib/3.8.2-gfbf-2023b
python3 tss_graph.py
