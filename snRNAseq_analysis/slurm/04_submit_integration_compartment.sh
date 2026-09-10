#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

for subtype in ER MED IMPC TNBC CYS APO PLC MpBC; do
  for compartment in Immune Stroma Tumour; do
    sbatch -J "${subtype}-${compartment}" \
      slurm/04_run_integration_compartment.sh "$subtype" "$compartment"
  done
done
