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

BAM_ROOT="/users/rg/baygun/BlaER_crg/files/long_read_files"
REF_DIR="/users/rg/baygun/BlaER_crg/long_read/epi2me/refdata"
RESULTS_BASE="/users/rg/baygun/BlaER_crg/long_read/epi2me/results"

mkdir -p "$RESULTS_BASE"

for RUN in t0 t120; do
  nextflow run epi2me-labs/wf-single-cell -r v3.3.0 \
    --expected_cells 10000 \
    --bam "${BAM_ROOT}/${RUN}_raw_files" \
    --kit '3prime:v3' \
    --ref_genome_dir "${REF_DIR}" \
    --out_dir "${RESULTS_BASE}/${RUN}_epi2me_result" \
    -profile singularity \
    -without-docker
done


conda deactivate
