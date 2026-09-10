#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

for subtype in ER MED IMPC TNBC CYS APO PLC MpBC; do
  for subpop in Myeloid T B_Plasma Fibroblast; do
    sbatch -J "${subtype}-${subpop}" \
      slurm/05_run_integration_subpopulation.sh "$subtype" "$subpop"
  done
done
