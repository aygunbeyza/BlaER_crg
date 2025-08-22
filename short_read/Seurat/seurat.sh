#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=seurat
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

module purge
module load R/4.4.2-gfbf-2023b

Rscript run_df_batch.R

