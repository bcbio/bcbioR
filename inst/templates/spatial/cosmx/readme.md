# Guidelines for analysis

Make sure there is a valid project name, and modify `information.R` with the right information for your project. You can use this file with any other Rmd to include the project/analysis information.

## QC

`QC/QC.Rmd` is a template for QC metrics. It plots the locations of cells on the slide, filters cells using the number of genes and AtoMX quality flags, and normalizes the data. It also provides sample code for clustering and cell type identification.  
 
Read instruction in the R and Rmd scripts to render it. 

<b>Note</b> that future versions of this template will include code for building your own RDS object out of .csv.gz files produced by the AtoMX software instead of loading an RDS directly, as this allows the analyst to access information about transcript locations and cell sementation that is not available in the pre-made RDS objects from AtoMx or BWH.

## DropBox

- In `reports/QC`
  - [ ] copy QC `Rmd/R/html/figures`
