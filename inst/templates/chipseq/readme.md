# Guidelines for analysis

Make sure there is a valid project name, and modify `information.R` with the right information for your project. You can use this file with any other Rmd to include the project/analysis information.

## Run data with nf-core rnaseq

This templates assume data has been processed by [nf-core/chipseq](https://nf-co.re/chipseq/2.1.0/docs/usage/).
We recommend to use the samplesheet.csv used with nf-core as metadata file, where other relevant columns can be there even if they are not used by the pipeline.

## QC

`QC/QC.Rmd` is a template for QC metrics. It includes basic read-level statistics, peak quality information, sample correlation analysis, and PCA that it produces using the above samplesheet and output from the nf-core pipeline. Use `params_qc.R` to provide the required input files. 

## DiffBind

`diffbind/diffbind.Rmd` is a template for comparing peak binding betweeen two groups. Use `params_diffbind.R` to provide the required input files. 

On the YAML header file of the Rmd you can specify some parameters including the conditions to be compared, the genome used, and the desired output file names. This template has examples of:
* calculating a peak counts matrix
* PCA
* differential binding analysis
* peak annotation
* functional analysis (coming soon)

This template writes to CSV a log2 normalized counts matrix of peaks x samples as well as the annotated significant results of the differential binding analysis. 
