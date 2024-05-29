# Guideline for RNAseq downstream analysis


## Downstream analysis

Please, modify `information.R` with the right information. You can use this file with any other Rmd to include the project/analysis information.

### QC

`QC/QC.Rmd` is a template for QC metrics. Use `params_qc.R` for `bcbio` 
 or `QC/QC_nf-core.Rmd` `params_qc_nf-core.R` for `nf-core/rnaseq` outputs.
 
Read instruction in the R and Rmd scripts to render it.

### DE

`DE/DEG.Rmd` is a template for two groups comparison. `params_de.R` has the information of the input files to load. You can point to `bcbio` or `nf-core/rnaseq` output files.

On the `YAML` header file of the `Rmd` you can specify some parameters or just set them up in the first chunk of code of the template. This template has examples of:

- sub-setting data
- two groups comparison
- volcano plot
- MA plot
- Pathway analysis
- Tables

There are some code related to alternative analysis:

- `DE/PCA_variance_analysis.R` that shows how to compare variance among groups to decide how to perform DE analysis.
- `DE/Multiplicative_DE_docs.md` that shows some cases when there is multiple variables in the model with multiple levels: sex (2 levels) and genotype (4 levels)
