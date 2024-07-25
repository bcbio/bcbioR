# Project name

## Start with cell-ranger

`pre-process-w-cellranger.md` contains step by step guidelines on how to run cellranger in O2 and load data into R. This `scripts/seurat_init.R` script contains all the pieces to go from cellranger output to Seurat obj. It is assuming a mouse genome.

## DropBox

-   In `reports/QC`
    -   [ ] copy QC `Rmd/R/html/figures`
-   In `reports/Clusters`
    -   [ ] the analysis of `SCTransform`, ,`RunPCA` ,`FindNeighbors`, ,`FindClusters`, `RunUMAP`
    - [ ] the analysis of `FindMarkers` and `Cell Identification`
-   In `reports/DE`, for *each analysis*:
    -   TBD

## GitHub

-   [ ] Push all `*Rmd` `*R` files used for the *QC* and *DE* analysis respecting folder structure.

Please, ignore `*html/figures/csv` and any output of the code.
