library(tidyverse)
library(SummarizedExperiment)
library(janitor)
load_metrics <- function(multiqc_data_dir){
  
  # the reading-in of these next two files needs changing in order to correctly 
  # account for samples that have been sequenced multiple times.
  # simply removing T1 is not the correct way to do it 
  fastqc <- read_tsv(file.path(multiqc_data_dir, 'multiqc_fastqc.txt')) %>% clean_names() %>%
    dplyr::select(sample, total_reads = total_sequences) %>%
    mutate(sample = gsub('_T1', '', sample))
  samtools <- read_tsv(file.path(multiqc_data_dir, 'multiqc_samtools_stats.txt')) %>% clean_names() %>%
    dplyr::select(sample, mapped_reads = reads_mapped) %>%
    mutate(sample = gsub('_T1', '', sample))
  
  
  phantom <- read_tsv(file.path(multiqc_data_dir, 'multiqc_phantompeakqualtools.txt')) %>% clean_names() %>%
    dplyr::select(sample, nsc, rsc)
  frip <- read_tsv(file.path(multiqc_data_dir, 'multiqc_frip_score-plot.txt')) %>% select(-Sample) %>% 
    pivot_longer(everything(), names_to = 'sample', values_to = 'frip') %>% filter(!is.na(frip))
  peak_count <- read_tsv(file.path(multiqc_data_dir, 'multiqc_peak_count-plot.txt')) %>% select(-Sample) %>% 
    pivot_longer(everything(), names_to = 'sample', values_to = 'peak_count') %>% filter(!is.na(peak_count))
  nrf <- read_tsv(file.path(multiqc_data_dir, 'mqc_picard_deduplication_1.txt')) %>% clean_names() %>%
    mutate(nrf = unique_unpaired / (unique_unpaired + duplicate_unpaired)) %>%
    dplyr::select(sample, nrf)
  
  metrics <- full_join(fastqc, samtools) %>% full_join(phantom) %>% full_join(frip) %>% 
    full_join(peak_count) %>% full_join(nrf) %>%
    mutate(mapped_reads_pct = round(mapped_reads/total_reads*100,1))
  
  metrics$sample <- make.names(metrics$sample)
  rownames(metrics) <- metrics$sample
  return(metrics)
}

load_coldata <- function(coldata_fn, column=NULL, numerator=NULL, denominator=NULL, subset_column = NULL, subset_value = NULL){
  coldata=read.csv(coldata_fn) %>%
    dplyr::distinct(sample, .keep_all = T) %>%
    dplyr::select(!matches("fastq")) %>%
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
  coldata$antibody <- ifelse(coldata$antibody == '', 'input', coldata$antibody)
  coldata$type <- ifelse(coldata$antibody == 'input', 'input', 'chip')

  if (!is.null(denominator))
    coldata[[column]] = relevel(as.factor(coldata[[column]]), denominator)

  return(coldata)
}

load_counts <- function(counts_fn){

  counts <- readRDS(counts_fn)
  return(counts)

}

load_peaks <- function(peaks_dir){
  peaks_fns <- list.files(peaks_dir, pattern = '_peaks.broadPeak')
  names(peaks_fns) <- gsub('_peaks.broadPeak', '', peaks_fns)
  peaks_all <- lapply(peaks_fns, function(fn) {
    peaks <- read_delim(file.path(peaks_dir, fn), col_names = F)
    peaks_df <- data.frame(peak_enrichment = peaks$X7, peak_rank = rank(dplyr::desc(peaks$X7))) %>% 
      dplyr::arrange(peak_rank)
    return(peaks_df)
  }) %>% bind_rows(.id = 'sample')
  return(peaks_all)
}
