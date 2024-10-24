# Guideline for RNAseq downstream analysis

Make sure there is a project name for this.

## Run data with nf-core rnaseq

This templates assume data has been processed by [nf-core/rnaseq](https://nf-co.re/rnaseq/3.14.0/docs/usage).
We recommend to use the samplesheet.csv used with nf-core as metadata file, where other relevant columns can be there even if they are not used by the pipeline.

## Downstream analysis

Modify `information.R` with the right information. You can use this file with any other Rmd to include the project/analysis information.

### QC

`QC/QC.Rmd` is a template for QC metrics. Use `params_qc.R` for `bcbio` 
 or `QC/QC_nf-core.Rmd` `params_qc_nf-core.R` for `nf-core/rnaseq` outputs.
 
Read instruction in the R and Rmd scripts to render it.

### DE

`DE/DEG.Rmd` is a template for comparison between two groups. `params_de.R` has the information for the input files to load. You can point to `bcbio` or `nf-core/rnaseq` output files.

On the `YAML` header file of the `Rmd` you can specify some parameters or just set them up in the first chunk of code of the template. This template has examples of:

- sub-setting data
- two groups comparison
- volcano plot
- MA plot
- Pathway analysis
- Tables

### Other templates

- `DE/GSVA.Rmd` shows an example on how to use [GSVA package](https://bioconductor.org/packages/release/bioc/html/GSVA.html) for estimating variation of gene set enrichment through the samples of a expression data set
- `DE/Cross-comparison-analysis.Rmd` shows an exmaple on how to compare two differential expression analysis from the `DEG.Rmd` template.
- `DE/Intersections.Rmd` shows an example on how to compare multiple differential expression analyses from `DE/DEG.Rmd` and find intersections.


