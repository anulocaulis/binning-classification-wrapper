#!/usr/bin/env bash
#SBATCH --job-name=vamb_assembly.megahit
#SBATCH --partition=math-alderaan
#SBATCH --account=biology-miller-annotation
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=6-23:30:00
#SBATCH --output=slurm_logs/manual_vamb_assembly.megahit_%j.out
#SBATCH --error=slurm_logs/manual_vamb_assembly.megahit_%j.err

set -euo pipefail
cd /storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper
bash scripts/manual_vamb/run_vamb_single_assembly.sh assembly.megahit
