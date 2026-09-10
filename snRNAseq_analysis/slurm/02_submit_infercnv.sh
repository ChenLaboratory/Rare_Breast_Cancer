#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

cut -f1 metadata/targets_all.txt | tail -n +2 | tr -d $'\r' | while read -r sample; do
  [ -z "$sample" ] && continue
  sbatch -J "${sample}_inf" slurm/02_run_infercnv.sh "$sample"
done
