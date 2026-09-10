#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=<your email>

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

module load mariadb-connector-c/3.3.10
module load geos/3.12.1
module load gdal/3.9.0
module load proj/9.4.0
module load hdf5/1.12.3
module load R/4.5.2

Rscript 05_Integration_Cell_CAF.R "$1"
