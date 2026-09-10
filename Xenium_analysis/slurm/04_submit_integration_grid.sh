#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

stage="${1:-main}"

for subtype in IMPC PLC CYS MED APO ER TNBC; do
  sbatch -J "$subtype" slurm/04_run_integration_grid.sh "$subtype" "$stage"
done
