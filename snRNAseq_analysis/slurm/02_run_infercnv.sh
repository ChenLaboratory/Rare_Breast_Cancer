#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=<your email>

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

JAGS_HOME=${JAGS_HOME:-"$HOME/jags"}
export PKG_CONFIG_PATH="${JAGS_HOME}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="${JAGS_HOME}/lib:${LD_LIBRARY_PATH:-}"

source slurm/_modules.sh

Rscript 02_InferCNV.R "$1"
