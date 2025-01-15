library(tidyverse)
library(SummarizedExperiment)
library(janitor)
load_metrics <- function(se_object, multiqc_data_dir, gtf_fn, counts){

  # bcbio input
  if (!is.na(se_object)){

    se <- readRDS(se_object)
    metrics <- metadata(se)$metrics %>% as.data.frame()
    # left_join(coldata %>% rownames_to_column('sample')) %>% column_to_rownames('sample')
  } else { #nf-core input

    # Get metrics from nf-core into bcbio like table
    # many metrics are already in the Genereal Table of MultiQC, this reads the file
    metrics <- read_tsv(file.path(multiqc_data_dir, 'multiqc_general_stats.txt'))

    # we get some more metrics from Qualimap and rename columns
    metrics_qualimap <- read_tsv(file.path(multiqc_data_dir, 'mqc_qualimap_genomic_origin_1.txt'))
    metrics <- metrics %>% full_join(metrics_qualimap)
    metrics <- metrics %>%
      clean_names() %>%
      dplyr::rename_with(~gsub('.*mqc_generalstats_', '', .))

    # This uses the fastqc metrics to get total reads
    total_reads <- metrics %>%
      dplyr::filter(!is.na(fastqc_raw_total_sequences)) %>%
      remove_empty(which = 'cols') %>%
      dplyr::rename(single_sample = sample) %>%
      mutate(sample = gsub('_[12]+$', '', single_sample)) %>%
      group_by(sample) %>%
      summarize(total_reads = sum(fastqc_raw_total_sequences))

    # This renames to user-friendly names the metrics columns
    metrics <- metrics %>%
      dplyr::filter(is.na(fastqc_raw_total_sequences)) %>%
      remove_empty(which = 'cols') %>%
      full_join(total_reads) %>%
      mutate(mapped_reads = samtools_reads_mapped) %>%
      mutate(exonic_rate = exonic/(star_uniquely_mapped * 2)) %>%
      mutate(intronic_rate = intronic/(star_uniquely_mapped * 2)) %>%
      mutate(intergenic_rate = intergenic/(star_uniquely_mapped * 2)) %>%
      mutate(x5_3_bias = qualimap_5_3_bias)

    # Sometimes we don't have rRNA due to mismatch annotation, We skip this if is the case
    gtf <- NULL
    if (genome =="other"){
      gtf <- gtf_fn
    }else{
      if (genome == "hg38") {
        gtf <- "hg38.rna.gtf.gz"
      } else if (genome == "mm10") {
        gtf <- "mm10.rna.gtf.gz"
      } else if (genome == "mm39") {
        gtf <- "mm39.rna.gtf.gz"
      }
      gtf <- system.file("extdata", "annotation",
                         gtf,
                         package="bcbioR")
    }
    if (is.null(gtf)) {
      warning("No genome provided! Please add it at the top of this Rmd")
    }

    gtf=rtracklayer::import(gtf)

    one=grep("gene_type", colnames(as.data.frame(gtf)), value = TRUE)
    another=grep("gene_biotype", colnames(as.data.frame(gtf)), value = TRUE)
    biotype=NULL
    if(length(one)==1){
      biotype=one
    }else if(length(another)==1){
      biotype=another
    }else{
      warning("No gene biotype founded")
    }
    metrics$sample <- make.names(metrics$sample)
    if (!is.null(biotype)){
      annotation=as.data.frame(gtf) %>% .[,c("gene_id", biotype)]
      rRNA=grepl("rRNA|tRNA",annotation[[biotype]])
      genes=intersect(annotation[rRNA,"gene_id"],row.names(counts))
      ratio=data.frame(sample=colnames(counts),
                       r_and_t_rna_rate=colSums(counts[genes,])/colSums(counts))
      metrics = left_join(metrics, ratio, by="sample")
    }else{
      metrics[["r_and_t_rna_rate"]] <- NA
    }

    # if ("custom_content_biotype_counts_percent_r_rna" %in% colnames(metrics)){
    #   metrics <- mutate(metrics, r_rna_rate = custom_content_biotype_counts_percent_r_rna)
    # }else{
    #  metrics[["r_rna_rate"]] <- NA
    # }
    metrics=metrics[,c("sample","mapped_reads","exonic_rate","intronic_rate",
                       "total_reads",
                       "x5_3_bias", "r_and_t_rna_rate","intergenic_rate")]
  }
  rownames(metrics) <- metrics$sample
  return(metrics)
}

load_coldata <- function(coldata_fn, column=NULL, subset_column = NULL, subset_value = NULL){
  coldata=read.csv(coldata_fn) %>%
    dplyr::distinct(sample, .keep_all = T) %>%
    dplyr::select(!matches("fastq"), !matches("strandness")) %>%
    distinct()
  if('description' %in% names(coldata)){
    coldata$sample <- tolower(coldata$description)
  }
  coldata <- coldata %>% distinct(sample, .keep_all = T)
  if (!is.null(column))
    stopifnot(column %in% names(coldata))

  # use only some samples, by default use all
  if (!is.null(subset_column)){
    coldata <- coldata[coldata[[paste(subset_column)]] == subset_value, ]
  }
  #coldata <- coldata[coldata[[paste(column)]] %in% c(numerator, denominator), ]
  #browser()
  coldata$sample <- make.names(coldata$sample)
  rownames(coldata) <- coldata$sample
  coldata$description <- coldata$sample

  # if (!is.null(denominator))
  #   coldata[[column]] = relevel(as.factor(coldata[[column]]), denominator)

  return(coldata)
}

load_counts <- function(counts_fn){

  # bcbio input
  if(grepl('csv', counts_fn)){
    counts <- read_csv(counts_fn) %>%
      mutate(gene = str_replace(gene, pattern = "\\.[0-9]+$", "")) %>%
      column_to_rownames('gene')
    colnames(counts) = tolower(colnames(counts))
    return(counts)
  } else { # nf-core input
    counts <- read_tsv(counts_fn) %>% dplyr::select(-gene_name) %>%
      mutate(gene_id = str_replace(gene_id, pattern = "\\.[0-9]+$", "")) %>%
      column_to_rownames('gene_id') %>% round %>% as.matrix()
    counts=counts[rowSums(counts)!=0,]
    return(counts)
  }

}
