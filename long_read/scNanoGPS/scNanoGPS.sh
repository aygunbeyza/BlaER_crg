#!/bin/bash
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --time=60:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=scNnaoGPS
#SBATCH --mail-user=beyza.aygun@crg.eu
#SBATCH --mail-type=FAIL

conda init
source ~/.bashrc
conda activate scNanoGPS

python scNanoGPS/run_scNanoGPS.py -i /T000.fastq.gz -d T0_scNanoGPS_res -p 3p -t 8 \
--gtf /annotation.gtf.gz \
--ref_genome /hg38.sorted.fa \
--idx_genome /hg38.mmi \
--tmp_dir T0_tmp

conda deactivate
