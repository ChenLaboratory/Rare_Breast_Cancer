#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")/.."

for subtype in ER MED IMPC TNBC CYS APO PLC MpBC; do
  sbatch -J "${subtype}-Final" slurm/06_run_integration_final.sh "$subtype"
done
