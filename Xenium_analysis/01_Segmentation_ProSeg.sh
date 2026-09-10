#!/bin/bash
##
## 01_Segmentation_ProSeg.sh

set -euo pipefail

sample=$1
run=$2
region=$3

## ---- Paths ----
XENIUM_DIR=${XENIUM_DIR:-"Xenium"}
PROSEG_DIR=${PROSEG_DIR:-"Proseg"}
SEG_DIR=${SEG_DIR:-"Proseg_mask"}
PROSEG_TO_BAYSOR=${PROSEG_TO_BAYSOR:-"proseg-to-baysor"}
XENIUMRANGER=${XENIUMRANGER:-"xeniumranger"}

mkdir -p "${SEG_DIR}/${sample}"
cd "${SEG_DIR}/${sample}"

## ---- ProSeg to Baysor ----
"${PROSEG_TO_BAYSOR}" \
    "${PROSEG_DIR}/${sample}/transcript-metadata.csv.gz" \
    "${PROSEG_DIR}/${sample}/cell-polygons.geojson.gz" \
    --output-transcript-metadata "${PROSEG_DIR}/${sample}/baysor-transcript-metadata.csv" \
    --output-cell-polygons "${PROSEG_DIR}/${sample}/baysor-cell-polygons.geojson"

## ---- Import the segmentation into the Xenium bundle ----
"${XENIUMRANGER}" import-segmentation \
    --id "${sample}" \
    --xenium-bundle "${XENIUM_DIR}/${run}/${region}" \
    --viz-polygons "${PROSEG_DIR}/${sample}/baysor-cell-polygons.geojson" \
    --transcript-assignment "${PROSEG_DIR}/${sample}/baysor-transcript-metadata.csv" \
    --units microns
