# Guideline for scRNAseq analysis

Make sure there is a valid project name, and modify `information.R` with the right information for your project. You can use this file with any other Rmd to include the project/analysis information.

# cell-ranger

`pre-process-w-cellranger.md` contains step by step guidelines on how to run cellranger and load data into R. This `scripts/seurat_init.R` script contains all the pieces to go from cellranger output to Seurat obj. It is assuming a mouse genome.

# QC

Currently we are working on deploying a shiny app to inspect the single cell object and find the best cut-offs for filtering. The Rmd that helps to visualize the before and after is `QC.Rmd`.

# Integration

`Integration/norm_integration.rmd` is a template with guidelines on how to work with multiple samples. It compares log2norm vs SCT, work with SCT by samples to remove batch biases better, provide options for integration between CCA and Harmony. As last step, it contains cell type clustering and visualization to help decide the best parameters.

# differential_expression

`scRNA_pseudobulk.Rmd` is a template that performs pseudobulk differential expression analysis using DESeq2. It takes as input:

-   information.R, containing basic facts about the project, analysis, and scientific question
-   Seurat object containing both raw counts and log normalized data
-   a metadata variable of interest and two levels of that variable to be compared
-   a column in the seurat object metadata containing cluster names/numbers
-   the name/number of the cluster of interest

The template aggregates single cell counts to the sample x cluster x factor_of_interest pseudobulk level, subsets to the cluster of interest, performs the desired statistical comparisons, and visualizes the results.
