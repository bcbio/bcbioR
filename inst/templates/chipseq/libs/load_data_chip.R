library(tidyverse)
library(SummarizedExperiment)
library(janitor)

load_metrics <- function(multiqc_data_dir){
  
  fastqc <- read_tsv(file.path(multiqc_data_dir, 'multiqc_fastqc.txt')) %>% clean_names() %>%
    dplyr::select(sample, total_reads = total_sequences) %>%
    mutate(new_sample = gsub('_T[0-9]+', '', sample)) %>%
    group_by(new_sample) %>% 
    summarize(new_total_reads = sum(total_reads)) %>%
    dplyr::select(sample = new_sample, total_reads = new_total_reads)
  samtools <- read_tsv(file.path(multiqc_data_dir, 'multiqc_samtools_stats.txt')) %>% clean_names() %>%
    dplyr::select(sample, mapped_reads = reads_mapped) %>%
    mutate(new_sample = gsub('_T[0-9]+', '', sample)) %>%
    group_by(new_sample) %>% 
    summarize(new_mapped_reads = sum(mapped_reads)) %>%
    dplyr::select(sample = new_sample, mapped_reads = new_mapped_reads)
  
  phantom <- read_tsv(file.path(multiqc_data_dir, 'multiqc_phantompeakqualtools.txt')) %>% clean_names() %>%
    dplyr::select(sample, nsc, rsc)
  frip <- read_tsv(file.path(multiqc_data_dir, 'multiqc_frip_score-plot.txt')) %>% dplyr::select(-Sample) %>% 
    pivot_longer(everything(), names_to = 'sample', values_to = 'frip') %>% filter(!is.na(frip))
  peak_count <- read_tsv(file.path(multiqc_data_dir, 'multiqc_peak_count-plot.txt')) %>% dplyr::select(-Sample) %>% 
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
  coldata <- read.csv(coldata_fn) %>%
    dplyr::distinct(sample_name, .keep_all = T) %>%
    dplyr::select(!matches("fastq")) %>%
    distinct()

  if (!is.null(column))
    stopifnot(column %in% names(coldata))

  # use only some samples, by default use all
  if (!is.null(subset_column)){
    coldata <- coldata[coldata[[paste(subset_column)]] == subset_value, ]
  }
 
  if (!is.null(denominator))
    coldata[[column]] = relevel(as.factor(coldata[[column]]), denominator)

  return(coldata)
}

load_counts <- function(counts_fn){

  counts <- readRDS(counts_fn)
  return(counts)

}

load_peaks <- function(peaks_dir){
  if(grepl('broadPeak', peaks_dir)){
    peaks_fns <- list.files(peaks_dir, pattern = '_peaks.broadPeak')
    names(peaks_fns) <- gsub('_peaks.broadPeak', '', peaks_fns)
  } else {
    peaks_fns <- list.files(peaks_dir, pattern = '_peaks.narrowPeak')
    names(peaks_fns) <- gsub('_peaks.narrowPeak', '', peaks_fns)
  }
  peaks_all <- lapply(peaks_fns, function(fn) {
    peaks <- read_delim(file.path(peaks_dir, fn), col_names = F)
    peaks_df <- data.frame(seqnames = peaks$X1, start = peaks$X2, end = peaks$X3,
                           peak_enrichment = peaks$X7, peak_rank = rank(dplyr::desc(peaks$X7))) %>% 
      dplyr::arrange(peak_rank)
    return(peaks_df)
  }) %>% bind_rows(.id = 'sample')
  peaks_all$sample_group <- gsub('_REP[0-9]+', '', peaks_all$sample)
  
  return(peaks_all)
}

make_diffbind_samplesheet <- function(coldata, bam_dir, peaks_dir, column = NULL){
  bam_files <- data.frame(bam = list.files(bam_dir, pattern = '.bam$', full.names = T)) %>%
    mutate(sample = sub("\\..*", "",basename(bam)))
  
  peak_files <- data.frame(Peaks = list.files(peaks_dir, pattern = 'Peak$', full.names = T)) %>%
    mutate(SampleID = sub("\\..*", "",basename(Peaks))) %>%
    mutate(SampleID = gsub('_peaks', '', SampleID))
  
  coldata_for_diffbind <- coldata %>% 
    dplyr::mutate(SampleID = paste0(sample, "_REP", replicate))

  coldata_for_diffbind$Condition <- coldata_for_diffbind[[column]]
  
  samplesheet <- coldata_for_diffbind %>% dplyr::select(SampleID, Condition, replicate) %>% 
    left_join(bam_files %>% dplyr::select(SampleID = sample, bamReads = bam), by = 'SampleID') %>%
    left_join(peak_files, by = 'SampleID')
  
  samplesheet$PeakCaller <- "narrow"
  
  return(samplesheet)
}

get_databases=function(sps="human"){
  all_in_life=list(
    msigdbr(species = sps, category = "H") %>% mutate(gs_subcat="Hallmark"),
    #  msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME"),
    msigdbr(species = sps, category = "C2", subcategory = "CP:KEGG"),
    #  msigdbr(species = "human", category = "C2", subcategory = "CP:PID"),
    msigdbr(species = sps, category = "C5", subcategory = "GO:BP"),
    msigdbr(species = sps, category = "C5", subcategory = "GO:MF")
    #  msigdbr(species = "human", category = "C5", subcategory = "HPO"),
    #  msigdbr(species = "human", category = "C3", subcategory = "TFT:GTRD"),
    #  msigdbr(species = "human", category = "C6") %>% mutate(gs_subcat="Oncogenic")
  )
  all_in_life
}

run_fora=function(input, uni,all_in_life){
  # browser()
  total_deg=length(unique(input))/length(unique(uni$ENTREZID))
  pathways_ora_all = lapply(all_in_life, function(p){
    pathway = split(x = p$entrez_gene, f = p$gs_name)
    db_name = paste(p$gs_cat[1], p$gs_subcat[1],sep=":")
    respath <- fora(pathways = pathway,
                    genes = unique(input$ENTREZID),
                    universe = unique(uni$ENTREZID),
                    minSize  = 15,
                    maxSize  = 500)
    # coll_respath = collapsePathwaysORA(respath[order(pval)][padj < 0.1],
    #                                    pathway, unique(input$ENTREZID), unique(uni$ENTREZID))
    as_tibble(respath)  %>%
      mutate(database=db_name, NES=(overlap/size)/(total_deg))
  }) %>% bind_rows() %>%
    mutate(analysis="ORA")
  ora_tb = pathways_ora_all %>% unnest(overlapGenes) %>%
    group_by(pathway) %>%
    left_join(uni, by =c("overlapGenes"="ENTREZID")) %>%
    dplyr::select(pathway, padj, NES, SYMBOL, analysis,
                  database) %>%
    group_by(pathway,padj,NES,database,analysis) %>%
    summarise(genes=paste(SYMBOL,collapse = ","))
  ora_tb
  
}

