# snRNA-seq Analysis Pipeline

Analysis workflow for the single nucleus RNA-seq dataset: per sample quality
control, copy number inference to separate tumour from normal nuclei, then a
three level integration of the samples of each subtype (all cells, then the
immune, stroma and tumour compartments, then the sub-populations of the immune
and stroma compartments), a final annotation that removes the contaminating
clusters, and pseudo-bulk differential expression between the subtypes.

Samples are grouped by subtype (APO, CYS, ER, IMPC, MED, MpBC, PLC, TNBC) in
`metadata/targets_<Subtype>.txt`.

## Workflow

Scripts, run in order:

| Script | Description |
|---|---|
| `00_Install_Packages.R` | R package installation |
| `01_Preprocessing_QC.R` | Per sample QC, doublet removal, clustering |
| `02_InferCNV.R` | Copy number inference against a normal reference, and the per nucleus instability score |
| `03_Integration_All.R` | All cell Harmony integration of one subtype, level 1 cell types |
| `04_Integration_Compartment.R` | Re-integration of one compartment (Immune, Stroma, Tumour), level 2 cell types |
| `05_Integration_Subpopulation.R` | Re-integration of one sub-population (Myeloid, T, B_Plasma, Fibroblast), level 3 clusters, GO and KEGG |
| `06_Integration_Final.R` | The three levels combined, tumour populations relabelled, contamination removed |
| `07_Pseudobulk_Tumour.R` | Pseudo-bulk DE and pathway analysis between the tumour cells of the subtypes |
| `08_Pseudobulk_CAF.R` | Pseudo-bulk DE, pathway analysis and signature scoring of the CAFs |

* `setup.R` defines the packages, the colour palettes, the cell type colour
  map used across every figure, the cell type groups that define the analysis
  hierarchy, and the EMT gene set;
* `config.R` defines the paths, the sample metadata readers, the per subtype
  clustering resolutions, the manual relabelling rules, and the helpers shared
  by more than one script.


## Directory structure

```text
<repository>/
|-- 00_Install_Packages.R
|-- setup.R
|-- config.R
|-- 01_Preprocessing_QC.R
|-- 02_InferCNV.R
|-- 03_Integration_All.R
|-- 04_Integration_Compartment.R
|-- 05_Integration_Subpopulation.R
|-- 06_Integration_Final.R
|-- 07_Pseudobulk_Tumour.R
|-- 08_Pseudobulk_CAF.R
|-- slurm/
|   |-- _modules.sh
|   |-- 01_submit_preprocessing.sh          01_run_preprocessing.sh
|   |-- 02_submit_infercnv.sh               02_run_infercnv.sh
|   |-- 03_submit_integration_all.sh        03_run_integration_all.sh
|   |-- 04_submit_integration_compartment.sh 04_run_integration_compartment.sh
|   |-- 05_submit_integration_subpopulation.sh 05_run_integration_subpopulation.sh
|   |-- 06_submit_integration_final.sh      06_run_integration_final.sh
|   |-- 07_run_pseudobulk_tumour.sh
|   |-- 08_run_pseudobulk_caf.sh
|-- metadata/
|   |-- targets_all.txt
|   |-- targets_<Subtype>.txt
|   |-- <Subtype>_celltype_{All,Immune,Stroma,Tumour}.txt
|   |-- Contamination.txt
|   |-- HsMarkers_snRNAseq.txt
|   |-- CAF_Signatures.txt
```


## Metadata

| File | Contents |
|---|---|
| `targets_all.txt` | One row per sample: the sequencing run (`Run`), the mitochondrial percentage cutoff (`MT`) and the clustering resolution of the per sample analysis (`Resolution`) |
| `targets_<Subtype>.txt` | Samples of one subtype, in the order used by the multi panel figures. One column, `Sample` |
| `<Subtype>_celltype_All.txt` | Manual cluster to level 1 cell type assignment |
| `<Subtype>_celltype_{Immune,Stroma,Tumour}.txt` | Manual cluster to level 2 cell type assignment, per compartment |
| `Contamination.txt` | Clusters judged to be contamination: `Subtype`, `Step` (the level they were seen at) and `Cluster` |
| `HsMarkers_snRNAseq.txt` | Marker gene panel, columns `Cell_Type` and `Genes` |
| `gene_ordering_table.tsv` | Gene positions required by inferCNV: symbol, chromosome, start, end |
| `CAF_Signatures.txt` | CAF signature gene lists, columns `CellType` and `Gene` |

The per subtype clustering resolutions, the manual tumour relabelling, the
contamination matching rules and the pseudo-bulk thresholds are all set in
`config.R`.

## Running the pipeline

On a SLURM cluster, from the repository root.

```bash
./slurm/01_submit_preprocessing.sh
./slurm/02_submit_infercnv.sh
./slurm/03_submit_integration_all.sh
./slurm/04_submit_integration_compartment.sh
./slurm/05_submit_integration_subpopulation.sh
./slurm/06_submit_integration_final.sh
sbatch slurm/07_run_pseudobulk_tumour.sh
sbatch slurm/08_run_pseudobulk_caf.sh
```

The same steps run one sample, subtype or compartment at a time:

```bash
Rscript 01_Preprocessing_QC.R SK01
Rscript 02_InferCNV.R SK01
Rscript 03_Integration_All.R ER
Rscript 04_Integration_Compartment.R ER Immune
Rscript 05_Integration_Subpopulation.R ER Myeloid
Rscript 06_Integration_Final.R ER
Rscript 07_Pseudobulk_Tumour.R
Rscript 08_Pseudobulk_CAF.R
```

