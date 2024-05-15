# Tipical steps for scRNAseq downstream analysis

# Set up

This template assumes there is 4 folders with FASTQ files.

- GEX
- ATAC
- ADT
- HTO

**NOTE** use absolute paths inside all these files to avoid errors.

## Check file names

FASTQ files should all start with the same sufix, cellrange will identify that as FASTQ from same samples. You may not have that, and you need to fix the filenames to looks like:

`20240417_GEX1_6_AT12013_S41_L001_R2_001.fastq.gz` -> `20240417_GEX1_AT12013_S41_L001_R2_001.fastq.gz`

The R scripts `scripts/fix_filenames.R` has some code to do that.

## Analysis of scRNAseq and scATACseq

`cellrange-arc` is used for this analysis. It needs a `libraries.csv` file. In the `meta` folder you will find a template. It needs the FASTQ_PATH to scRNAseq and scATACseq with filenames corrected.

`scripts/gex_atac.sbatch` is an example of how to run this. It is assuming using HMS O2 cluster.

## Analysis of scADT with Hashing

`cellranger multi` is used for this analysis. It needs two config files:

- `feature_reference.csv` is a list of barcodes used to identify the proteins and the hashing. The example file in `meta` folder shows how to add 7 different hashing barcodes. 
- `multiconfig.csv` has three sections, to point to reference genome, input files (scRNAseq, scADT, HTO), the `feature_reference.csv` and samples. The example file in `meta` folder shows an example with 7 hasjing barcodes.

`scripts/gex_adt_hto.sbatch` is an example of how to run this. It is assuming using HMS O2 cluster.

# DropBox

-   In `reports/QC`
    -   [ ] copy QC `Rmd/R/html/figures`
-   In `reports/Clusters`
    -   [ ] the analysis of `SCTransform`, ,`RunPCA` ,`FindNeighbors`, ,`FindClusters`, `RunUMAP`
    - [ ] the analysis of `FindMarkers` and `Cell Identification`
-   In `reports/DE`, for *each analysis*:
    -   TBD

# GitHub

-   [ ] Push all `*Rmd` `*R` files used for the *QC* and *DE* analysis respecting folder structure.

Please, ignore `*html/figures/csv` and any output of the code.
