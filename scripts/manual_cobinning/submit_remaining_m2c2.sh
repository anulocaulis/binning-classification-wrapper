#!/usr/bin/env bash
set -euo pipefail

cd /storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper

sbatch scripts/manual_cobinning/idbaud_metabat2.sbatch.sh
sbatch scripts/manual_cobinning/idbaud_concoct.sbatch.sh
sbatch scripts/manual_cobinning/megahit_metabat2.sbatch.sh
sbatch scripts/manual_cobinning/megahit_concoct.sbatch.sh
sbatch scripts/manual_cobinning/metaspades_metabat2.sbatch.sh
sbatch scripts/manual_cobinning/metaspades_concoct.sbatch.sh
