# Guideline for scRNAseq analysis

Make sure there is a valid project name, and modify `information.R` with the right information for your project. You can use this file with any other Rmd to include the project/analysis information.

# cell-ranger

`pre-process-w-cellranger.md` contains step by step guidelines on how to run cellranger and load data into R. This `scripts/seurat_init.R` script contains all the pieces to go from cellranger output to Seurat obj. It is assuming a mouse genome.

# QC

Currently we are working on deploying a shiny app to inspect the single cell object and find the best cut-offs for filtering. The Rmd that helps to visualize the before and after is `QC.Rmd`.

# Integration

`Integration/norm_integration.rmd` is a template with guidelines on how to work with multiple samples. It compares log2norm vs SCT, work with SCT by samples to remove batch biases better, provide options for integration between CCA and Harmony. As last step, it contains cell type clustering and visualization to help decide the best parameters.

