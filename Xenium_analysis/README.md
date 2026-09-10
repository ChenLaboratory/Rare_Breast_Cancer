# Xenium Spatial Analysis Pipeline

Analysis workflow for the Xenium in situ dataset: cell segmentation with
ProSeg, per sample cell type annotation and hexagonal bin (niche) analysis,
grid level integration of the samples of each subtype into common niches, and
cell level integration of the fibroblast and CAF compartment.

Samples and parameters are listed in `metadata/targets_all.txt`.

## Workflow

Scripts, run in order:

| Script | Description |
|---|---|
| `00_Install_Packages.R` | R package installation |
| `01_Segmentation_ProSeg.sh` | Converts the ProSeg segmentation to Baysor format and imports it into the Xenium bundle |
| `02_Individual_Analysis.R` | Per sample cell clustering, manual cell type annotation, and niche analysis |
| `03_Spatial_LR.R` | Optional: spatial ligand receptor co-expression analysis |
| `04_Integration_Grid.R` | Grid level Harmony integration into common niches across the samples of a subtype |
| `05_Integration_Cell_CAF.R` | Cell level Harmony integration of the fibroblasts and CAFs of a subtype |

* `setup.R` defines the packages, the colour palettes, the cell type colour
  map used across every figure, and the marker gene sets;
* `config.R` defines the paths, the sample metadata readers, the per subtype
  analysis parameters, the sample specific special cases, and the helpers
  shared by more than one script.

Scripts 02 and 04 run twice. The `main` stage uses the manual cell type
annotation; the `post_CAF` stage repeats the same analysis after the CAF
clusters of script 05 have been merged back into the cell type labels, and
writes its objects to a separate directory.

## Directory structure

```text
<repository>/
|-- 00_Install_Packages.R
|-- setup.R
|-- config.R
|-- 01_Segmentation_ProSeg.sh
|-- 02_Individual_Analysis.R
|-- 03_Spatial_LR.R
|-- 04_Integration_Grid.R
|-- 05_Integration_Cell_CAF.R
|-- slurm/
|   |-- 01_submit_segmentation.sh        01_run_segmentation.sh
|   |-- 02_submit_individual.sh          02_run_individual.sh
|   |-- 03_submit_spatial_lr.sh          03_run_spatial_lr.sh
|   |-- 04_submit_integration_grid.sh    04_run_integration_grid.sh
|   |-- 05_submit_integration_cell_caf.sh 05_run_integration_cell_caf.sh
|-- metadata/
|   |-- targets_proseg_all.txt
|   |-- targets_all.txt
|   |-- targets_xenium_<Subtype>.txt
|   |-- <Subtype>_celltype_Individual.txt
|   |-- <Subtype>_caf_annotation.txt
|   |-- HsMarkers_Xenium.txt
```


## Metadata

| File | Contents |
|---|---|
| `targets_proseg_all.txt` | Sample, Run and Region of every Xenium output bundle |
| `targets_all.txt` | One row per sample: the ProSeg run, the cell filtering cutoffs (`hcut` genes, `vcut` transcripts), the cell and bin clustering resolutions (`resCell`, `resGrid`), the figure dimensions, the tissue splitting line (`Grad`, `Int`), and the minimum bin library size (`cutGrid`) |
| `targets_xenium_<Subtype>.txt` | Samples of one subtype, in the order used by the multi panel figures |
| `<Subtype>_celltype_Individual.txt` | Manual cluster to cell type assignment, one row per sample and cluster |
| `<Subtype>_caf_annotation.txt` | Manual CAF cluster to CAF subtype assignment, one row per CAF cluster of script 05 |
| `HsMarkers_Xenium.txt` | Marker gene panel used for the signature heatmap of script 04 |

Per subtype clustering resolutions, the artefactual niches excluded from the
integration, the pseudo-bulk library size cutoffs, and the sample specific
special cases are all set in `config.R`.

## Running the pipeline

On a SLURM cluster, from the repository root:

```bash
./slurm/01_submit_segmentation.sh
./slurm/02_submit_individual.sh
./slurm/04_submit_integration_grid.sh
./slurm/05_submit_integration_cell_caf.sh
```

Then repeat scripts 02 and 04 with the refined CAF labels:

```bash
./slurm/02_submit_individual.sh post_CAF
./slurm/04_submit_integration_grid.sh post_CAF
```

The optional ligand receptor analysis can be run at any point after script 02:

```bash
./slurm/03_submit_spatial_lr.sh
```

