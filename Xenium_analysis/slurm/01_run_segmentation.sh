#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=<your email>

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

./01_Segmentation_ProSeg.sh "$1" "$2" "$3"
