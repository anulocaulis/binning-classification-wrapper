#!/usr/bin/env bash
#SBATCH --job-name=m2c2_idbaud_metabat2
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=6-23:30:00
#SBATCH --output=slurm_logs/manual_idbaud_metabat2_%j.out
#SBATCH --error=slurm_logs/manual_idbaud_metabat2_%j.err

set -euo pipefail
cd /storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper
bash scripts/manual_cobinning/run_metawrap_single_binner.sh assembly.idbaud metabat2
