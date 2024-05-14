## code to prepare `DATASET` dataset goes here
library(DESeq2)
library(tidyverse)
coldata_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/coldata.csv'
counts_fn = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/salmon.merged.gene_counts.tsv'
multiqc_data_dir = 'https://raw.githubusercontent.com/bcbio/bcbioR-test-data/main/rnaseq/nf-core/multiqc-report-data/'
se_object = NA
column = "sample_type"
numerator = "tumor"
denominator = "normal"
subset_column = NULL
subset_value = NULL

coldata <- load_coldata(coldata_fn, column,
                        numerator, denominator,
                        subset_column, subset_value)
metrics <- load_metrics(se_object, multiqc_data_dir) %>%
  left_join(coldata, by = c('sample' = 'description')) %>%
  as.data.frame()
rownames(metrics) <- metrics$sample
counts <- load_counts(counts_fn)
counts <- counts[,colnames(counts) %in% coldata$description]
# if the names don't match in order or string check files names and coldata information
counts = counts[rownames(coldata)]
stopifnot(all(names(counts) == rownames(coldata)))

dds_to_use <- DESeqDataSetFromMatrix(counts, coldata, design = ~sample_type)
bcbio_vsd_data <- vst(dds_to_use)
usethis::use_data(bcbio_vsd_data, overwrite = TRUE)
