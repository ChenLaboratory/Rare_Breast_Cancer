#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

targets="metadata/targets_all.txt"

cut -f1 "$targets" | tail -n +2 | tr -d $'\r' | while read -r sample; do
  [ -z "$sample" ] && continue
  sbatch -J "LR_$sample" slurm/03_run_spatial_lr.sh "$sample"
done
