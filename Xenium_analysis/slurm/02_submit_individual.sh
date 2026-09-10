#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

stage="${1:-main}"
targets="metadata/targets_all.txt"

cut -f1 "$targets" | tail -n +2 | tr -d $'\r' | while read -r sample; do
  [ -z "$sample" ] && continue
  sbatch -J "$sample" slurm/02_run_individual.sh "$sample" "$stage"
done
