#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

for subtype in IMPC PLC CYS MED APO ER TNBC; do
  sbatch -J "${subtype}_CAF" slurm/05_run_integration_cell_caf.sh "$subtype"
done
