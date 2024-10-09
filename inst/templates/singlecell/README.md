# Project name

# Start with cell-ranger

`pre-process-w-cellranger.md` contains step by step guidelines on how to run cellranger and load data into R. This `scripts/seurat_init.R` script contains all the pieces to go from cellranger output to Seurat obj. It is assuming a mouse genome.

# QC

Currently we are working on deploying a shiny app to inspect the single cell object and find the best cut-offs for filtering. The Rmd that helps to visualize the before and after is `QC.Rmd`.

# Integration

Currently we are working on guidelines and templates for this step. There is some draft under *Integration** folder.

