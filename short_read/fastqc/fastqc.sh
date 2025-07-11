#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=fastqc
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

conda init
source ~/.bashrc
conda activate qc3

fastqc /users/rg/baygun/BlaER/files/short_read_files/*.fastq.gz -o /users/rg/baygun/BlaER/short_read/fastqc/result_fastqc

conda deactivate
