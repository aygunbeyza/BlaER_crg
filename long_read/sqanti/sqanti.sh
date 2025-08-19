#!/bin/bash
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sqanti
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

source ~/.bashrc
conda activate sqanti

SQANTI=/users/rg/baygun/BlaER_crg/long_read/sqanti/SQANTI3/sqanti3_qc.py
REFGTF=/users/rg/baygun/BlaER_crg/long_read/sqanti/ref/gencode.v47.primary_assembly.annotation.gtf
REFFA=/users/rg/baygun/BlaER_crg/long_read/sqanti/ref/hg38.sorted.fa
OUTROOT=/users/rg/baygun/BlaER_crg/long_read/sqanti/results

# lyric files
LYRIC_DIR=/users/rg/baygun/BlaER_crg/long_read/sqanti/ref
LYRIC_T0=$(ls -1 ${LYRIC_DIR}/*.gff | grep -E 'Hpre' | head -n1)
LYRIC_T120=$(ls -1 ${LYRIC_DIR}/*.gff | grep -E 'Hv3' | head -n1)

# isoquant files
ISO_T0="/users/rg/baygun/BlaER_crg/long_read/IsoQuant/result_isoquant_t0/OUT/OUT.transcript_models.gtf"
ISO_T120="/users/rg/baygun/BlaER_crg/long_read/IsoQuant/result_isoquant_t120/OUT/OUT.transcript_models.gtf"

declare -A RUNS=(
  [lyric_t0]="$LYRIC_T0"
  [lyric_t120]="$LYRIC_T120"
  [isoquant_t0]="$ISO_T0"
  [isoquant_t120]="$ISO_T120"
)

mkdir -p "$OUTROOT"

for NAME in "${!RUNS[@]}"; do
  ISOFORMS="${RUNS[$NAME]}"
  mkdir -p "$OUTROOT/$NAME"
  echo ">> Running SQANTI3 for $NAME"

  python "$SQANTI" \
    --isoforms "$ISOFORMS" \
    --refGTF "$REFGTF" \
    --refFasta "$REFFA" \
    -o "$NAME" -d "$OUTROOT/$NAME" \
    --cpus 4 --report both
done


conda deactivate
