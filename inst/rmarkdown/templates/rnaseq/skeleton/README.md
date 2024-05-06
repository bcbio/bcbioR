# Guideline for RNAseq downstream analysis

Make sure there is a project name for this.

## Run data with nf-core rnaseq

- Make sure you have access to our [Seqera WorkSpace](https://cloud.seqera.io/orgs/HBC/workspaces/core_production/launchpad)
- Transfer data to HCBC S3: Ask Alex/Lorena. Files will be at our S3 bucket `input/rawdata` folder
- Prepare the CSV file according this [instructions](https://nf-co.re/rnaseq/3.14.0/docs/usage#multiple-runs-of-the-same-sample). File should look like this:

```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,s3path/AEG588A1_S1_L002_R1_001.fastq.gz,s3path/AEG588A1_S1_L002_R2_001.fastq.gz,auto
CONTROL_REP1,s3path/AEG588A1_S1_L003_R1_001.fastq.gz,s3path/AEG588A1_S1_L003_R2_001.fastq.gz,auto
CONTROL_REP1,s3path/AEG588A1_S1_L004_R1_001.fastq.gz,s3path/AEG588A1_S1_L004_R2_001.fastq.gz,auto
```
- Upload file to our `Datasets` in Seqera using the name of the project but starting with `nfcore-rnaseq`
- Go to `Launchpad`, select `nf-core_rnaseq` pipeline, and select the previous created `Datasets` in the `input` parameter after clicking in `Browser`
  - Select an output directory with the same name used for the `Dataset` inside the `results` folder in S3
- When pipeline is down, data will be copied to our on-premise HPC

## Downstream analysis

## DropBox

- In `reports/QC`
  - [ ] copy `bcbio-se.rds` and `tximport-counts.csv`
  - [ ] copy QC `Rmd/R/html/figures`
- In `reports/DE`
  -	[ ] Normalized counts for all genes x all samples (csv format)
- In `reports/DE`, for *each analysis*:
  - **Note** For multiple comparisons/analysis, do a single report/template if possible in the parent folder using parameters whenever possible. 
  - Create a folder with the comparison names in the files. Numbering by comparison (`01.1_DE_comp1`, `01.2_DE_comp2`, etc.). If you’re running multiple models for the same comparison, append `_M#`. Add the following files under each folder:
  - [ ] Normalized count table with the samples used in this analysis/comparison.
  -	[ ] Full results `DESeq2` for all genes (csv format) with annotation columns appended. 
  -	[ ] Significant genes results file (subset of annotated full results by chosen p-value and LFC). Separate files will be created for each individual contrast.
  -	[ ] Significant genes results file as described above, but additionally append columns containing normalized count values for each sample.
  -	Make sure to append the gene symbols to these tables so the researcher can interpret the results.

## GitHub

- [ ] Push all `*Rmd` `*R` files used for the *QC* and *DE* analysis respecting folder structure.

Please, ignore `*html/figures/csv` and any output of the code.
