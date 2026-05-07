#!/usr/bin/env bash
set -euo pipefail

cd /storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper
mkdir -p slurm_logs

sbatch scripts/manual_vamb/flye_vamb.sbatch.sh
sbatch scripts/manual_vamb/idbaud_vamb.sbatch.sh
sbatch scripts/manual_vamb/megahit_vamb.sbatch.sh
sbatch scripts/manual_vamb/metaconnet_vamb.sbatch.sh
sbatch scripts/manual_vamb/metamdbg_vamb.sbatch.sh
sbatch scripts/manual_vamb/metaspades_vamb.sbatch.sh
sbatch scripts/manual_vamb/metaspades_hybrid_vamb.sbatch.sh
