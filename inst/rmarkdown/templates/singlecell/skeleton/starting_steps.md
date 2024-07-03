---
title: "From raw data to Seurat"
---


# Overview

This tutorial assumes that you are starting with 10x genomic data that has not yet been run through cellranger. If you have output files from cellranger (raw_feature_bc_matrix.h5) files skip to step 2. 

# Step 1 running cellranger

## Set up

Here are the steps that need to be completed prior to running cellranger.  

### Locate or create your genome

#### I have mouse or human

We have prebuilt references for mouse and human located here:

#### I have another genome

It is easy to generate a cellranger reference for any genome. All you need as input are a fasta file and a gtf file. Here is some [information](https://kb.10xgenomics.com/hc/en-us/articles/115003327112-How-can-we-add-genes-to-a-reference-package-for-Cell-Ranger) about what is required for the gtf file.
 
**Note: what is listed as "gene_id" (required in gtf) or "gene_name" (if used will be preferred) will be your row names (i.e. gene names). Make sure this is something useful or can be connected to information on what these genes are.**

Below is an example script for a non-model reference

```
#!/bin/sh
#SBATCH --partition=short
#SBATCH -o run.o
#SBATCH -e run.e
#SBATCH -t 0-2:30
#SBATCH -c 1
#SBATCH --mem=48G

module load cellranger/7.1.0

cellranger mkref \
    --genome=my_nonmodel_genome \
    --fasta=/path/to/my/genome/fasta/file/my_nonmodel_genome.fasta \
    --genes=/path/to/my/genome/annotation/file/my_nonmodel_genome.gtf

```


### Fastq data

You should have a number of output files from each sample. These should look like those below:

```
sample1_I1_001.fastq.gz 
sample1_I2_001.fastq.gz 
sample1_R1_001.fastq.gz 
sample1_R2_001.fastq.gz 
```

Cellranger will be looking for both the I and R files. If you do not have both you may have to run demultiplexing, [See here](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-mkfastq).

**NOTE: you may have multiple lanes per sample. there is no need to concatentate these prior to running cellranger.**


It is best to create one folder per sample with that sample name and put all of the files there.

**Cellranger expects 1 folder per sample**

Here is an example file structure

```
fastq_files
├── sample1
│   ├── sample1_I1_001.fastq.gz
│   ├── sample1_I2_001.fastq.gz
│   ├── sample1_R1_001.fastq.gz
│   ├── sample1_R2_001.fastq.gz
├── sample2
│   ├── sample2_I1_001.fastq.gz
```

## Run Cellranger

The easiest way to run cellranger is using the array feature on O2. [Here](https://github.com/hbc/knowledgebase/blob/master/rc/arrays_in_slurm.md) is a tutorial on arrays.

To run cellranger as an array you will need one extra file. This file called `samples.txt` will have the name of each sample on its own line.

```
sample1
sample2
sample3
...
sampleN
```

for ease `samples.txt` should be in the same directory as your sbatch script.

Here is an example sbatch script for running cellranger as an array

```(bash)
#!/bin/bash

#SBATCH --job-name=CellRangerCount3      # Job name
#SBATCH --partition=short            # Partition name
#SBATCH --time=0-05:00                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=16             # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G                     # Memory needed per node (total)
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID


samp=$(awk -v  awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' samples.txt)  ### This line will take the numeric slurm array task id and find the corresponding line number in samples.txt. The sample name is made into a variable called samp.

module load cellranger/7.1.0

cellranger count \
    --id=${samp} \ ## This is what your output folders will be named
    --fastqs=/path/to/your/folder/of/fastq/folders/${samp} \
    --transcriptome=/path/to/your/genome \
    --localcores=16 \
    --localmem=128 

```

This script can be run depending on the number of samples you have. Here we will call it N:

```
sbatch --array=1-N run_cellranger.sh
```

For example if you have 9 samples to run you can use:

```
sbatch --array=1-9 run_cellranger.sh
```

Arrays are also handy if you need to re-run just a single sample. Let's say you need to re-run the 1st and 9th sample in samples.txt

```
sbatch --array=1,9 run_cellranger.sh
```

# Step 2 - going from cellranger output to Seurat

