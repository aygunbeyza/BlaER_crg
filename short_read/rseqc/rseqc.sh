#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=rseqc
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

conda init
source ~/.bashrc
conda activate rseqc

module purge
module load R/4.4.2-gfbf-2023b

ref="annotation.bed12"
input_dir="/star"

# BAM files
bams=(
  "t0_Aligned.sortedByCoord.out.bam"
  "t0_with_genfull_Aligned.sortedByCoord.out.bam"
  "t120_Aligned.sortedByCoord.out.bam"
  "t120_with_genefull_Aligned.sortedByCoord.out.bam"
)

for bam in "${bams[@]}"; do
  junction_saturation.py \
    -i "$input_dir/$bam" \
    -r "$ref" \
    -o "${bam%.bam}"
done

conda deactivate
