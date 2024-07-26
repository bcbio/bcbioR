# Project name

# Start with cell-ranger

`pre-process-w-cellranger.md` contains step by step guidelines on how to run cellranger in O2 and load data into R. This `scripts/seurat_init.R` script contains all the pieces to go from cellranger output to Seurat obj. It is assuming a mouse genome.

# QC

Currently we are working on deploying a shiny app to inspect the single cell object and find the best cut-offs for filtering. The Rmd that helps to visualize the before and after is `QC.Rmd`.

# Integration

Currently we are working on guidelines and templates for this step. There is some draft under *Integration** folder.

# Cell to cell communication

CellChat template is at `CellToCell/cellchat.Rmd`. We have built a stable environment in O2 using the following modules:

```
# gcc/9.2.0 imageMagick/7.1.0 geos/3.10.2 cmake/3.22.2 R/4.3.1 fftw/3.3.10 gdal/3.1.4 udunits/2.2.28  boost/1.75.0 python/3.9.14
```
