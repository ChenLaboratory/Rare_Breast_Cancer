#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

for subtype in ER MED IMPC TNBC CYS APO PLC MpBC; do
  sbatch -J "${subtype}-All" slurm/03_run_integration_all.sh "$subtype"
done
