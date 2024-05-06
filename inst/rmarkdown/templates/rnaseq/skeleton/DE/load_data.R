library(tidyverse)
library(SummarizedExperiment)

load_metrics <- function(se_object, multiqc_data_dir){

  # bcbio input
  if (!is.na(se_object)){

    se <- readRDS(se_object)
    metrics <- metadata(se)$metrics %>% mutate(sample = toupper(sample)) %>% as.data.frame()
      # left_join(coldata %>% rownames_to_column('sample')) %>% column_to_rownames('sample')

  } else { #nf-core input

    metrics <- read_tsv(paste0(multiqc_data_dir, 'multiqc_general_stats.txt'))
    metrics_qualimap <- read_tsv(paste0(multiqc_data_dir, 'mqc_qualimap_genomic_origin_1.txt'))

    metrics <- metrics %>% full_join(metrics_qualimap)
    metrics <- metrics %>%
      clean_names() %>%
      dplyr::rename_with(~gsub('.*mqc_generalstats_', '', .))

    total_reads <- metrics %>%
      filter(grepl('_', sample)) %>%
      remove_empty(which = 'cols') %>%
      dplyr::rename(single_sample = sample) %>%
      mutate(sample = gsub('_[12]+$', '', single_sample)) %>%
      group_by(sample) %>%
      summarize(total_reads = sum(fastqc_raw_total_sequences))

    metrics <- metrics %>%
      filter(!grepl('_', sample)) %>%
      remove_empty(which = 'cols') %>%
      full_join(total_reads) %>%
      mutate(mapped_reads = samtools_reads_mapped) %>%
      mutate(exonic_rate = exonic/(star_uniquely_mapped * 2)) %>%
      mutate(intronic_rate = intronic/(star_uniquely_mapped * 2)) %>%
      mutate(intergenic_rate = intergenic/(star_uniquely_mapped * 2)) %>%
      mutate(r_rna_rate = custom_content_biotype_counts_percent_r_rna) %>%
      mutate(x5_3_bias = qualimap_5_3_bias) %>%
      dplyr::select(
        sample,
        total_reads,
        mapped_reads,
        exonic_rate,
        intronic_rate,
        intergenic_rate,
        r_rna_rate,
        x5_3_bias
      ) %>% as.data.frame()

    # metrics <- metrics %>%
    #   full_join(meta_df , by = c("sample" = "sample"))
  }
  rownames(metrics) <- metrics$sample
  return(metrics)
}

load_coldata <- function(coldata_fn, column, numerator, denominator, subset_column = NULL, subset_value = NULL){
  coldata=read.csv(coldata_fn)
  stopifnot(column %in% names(coldata))

  # use only some samples, by default use all
  if (!is.null(subset_column)){
    coldata <- coldata[coldata[[paste(subset_column)]] == subset_value, ]
  }
  #coldata <- coldata[coldata[[paste(column)]] %in% c(numerator, denominator), ]

  if('description' %in% names(coldata)){
    rownames(coldata) <- coldata$description
  }else {
    rownames(coldata) <- coldata$sample
    coldata$description <- coldata$sample
  }
  #browser()
  coldata[[column]] = relevel(as.factor(coldata[[column]]), denominator)

  return(coldata)
}

load_counts <- function(counts_fn){

  # bcbio input
  if(grepl('csv', counts_fn)){
    counts <- read_csv(counts_fn) %>% column_to_rownames('gene')
  } else { # nf-core input
    counts <- read_tsv(counts_fn) %>% dplyr::select(-gene_name) %>%
      mutate(gene_id = str_sub(gene_id, 1, 15)) %>%
      column_to_rownames('gene_id') %>% round

    return(counts)
  }

}
