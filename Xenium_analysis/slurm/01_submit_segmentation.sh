#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

targets="metadata/targets_proseg_all.txt"

tail -n +2 "$targets" | while IFS=$'\t' read -r sample run region; do
  region="${region//$'\r'/}"
  sbatch -J "$sample" slurm/01_run_segmentation.sh "$sample" "$run" "$region"
done
