#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=<your email>

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source slurm/_modules.sh

Rscript 05_Integration_Subpopulation.R "$1" "$2"
