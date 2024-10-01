# Guidelines for analysis

Make sure there is a valid project name, and modify `information.R` with the right information for your project. You can use this file with any other Rmd to include the project/analysis information.

## Run data with nf-core rnaseq

- Make sure you have access to our [Seqera WorkSpace](https://cloud.seqera.io/orgs/HBC/workspaces/core_production/launchpad)
- Transfer data to HCBC S3: Ask Alex/Lorena. Files will be at our S3 bucket `input` folder
- Prepare the CSV file according these [instructions](https://nf-co.re/chipseq/2.0.0/docs/usage/). File should look like this:

```csv
sample,fastq_1,fastq_2,antibody,control
WT_BCATENIN_IP_REP1,BLA203A1_S27_L006_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP_REP2,BLA203A25_S16_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP_REP2,BLA203A25_S16_L002_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP_REP2,BLA203A25_S16_L003_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_BCATENIN_IP_REP3,BLA203A49_S40_L001_R1_001.fastq.gz,,BCATENIN,WT_INPUT
WT_INPUT_REP1,BLA203A6_S32_L006_R1_001.fastq.gz,,,
WT_INPUT_REP2,BLA203A30_S21_L001_R1_001.fastq.gz,,,
WT_INPUT_REP2,BLA203A30_S21_L002_R1_001.fastq.gz,,,
WT_INPUT_REP3,BLA203A31_S21_L003_R1_001.fastq.gz,,,
```

You can add more columns to this file with more metadata, and use this file as the `coldata` file in the templates.

- Upload file to our `Datasets` in Seqera using the name of the project but starting with `chipseq-pi_lastname-hbc_code`
- Go to `Launchpad`, select `nf-core_chipseq` pipeline, and select the previous created `Datasets` in the `input` parameter after clicking in `Browser`
  - Select an output directory with the same name used for the `Dataset` inside the `results` folder in S3
- When pipeline is done, data will be copied to our on-premise HPC in the scratch system under `scratch/groups/hsph/hbc/bcbio/` folder

## QC

`QC/QC.Rmd` is a template for QC metrics. It includes basic read-level statistics, peak quality information, sample correlation analysis, and PCA that it produces using the above samplesheet and output from the nf-core pipeline. Use `params_qc.R` to provide the required input files. 

## DiffBind

`diffbind/diffbind.Rmd` is a template for comparing peak binding betweeen two groups. Use `params_diffbind.R` to provide the required input files. 

On the YAML header file of the Rmd you can specify some parameters including the conditions to be compared, the genome used, and the desired output file names. This template has examples of:
* calculating a peak counts matrix
* PCA
* differential binding analyiss
* peak annotation
* functional analysis (coming soon)

