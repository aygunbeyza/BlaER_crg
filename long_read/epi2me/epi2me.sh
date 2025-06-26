#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=64
#SBATCH --job-name=epi2me
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

conda init
source ~/.bashrc
conda activate qc3
module load Nextflow/24.10.3

nextflow run epi2me-labs/wf-single-cell -r v3.3.0 \
--expected_cells 10000 \
--bam /users/rg/baygun/BlaER_crg/files/long_read_files/t120_raw_files \
--kit '3prime:v3' \
--ref_genome_dir /users/rg/baygun/epi2me/refdata \
--out_dir /users/rg/baygun/epi2me/results_epi2me_120 \
-profile standard,singularity \
-without-docker

conda deactivate
