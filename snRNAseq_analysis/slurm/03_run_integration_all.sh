#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=<your email>

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source slurm/_modules.sh

Rscript 03_Integration_All.R "$1"
