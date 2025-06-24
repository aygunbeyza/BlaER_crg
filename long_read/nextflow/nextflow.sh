#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=nextflow_nanoplot_4
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL


conda init
source ~/.bashrc
conda activate qc3
module load Nextflow/24.10.3
nextflow run main.nf --outdir /users/rg/baygun/BlaER_crg/long_read/nextflow/result_nanoplot

conda deactivate
