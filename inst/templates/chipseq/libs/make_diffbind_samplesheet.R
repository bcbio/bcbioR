library(tidyverse)
library(tools)
library(DiffBind)
library(qs)
source('../libs/load_data.R')
source('params_diffbind-example.R')

coldata <- load_coldata(coldata_fn)

bam_files <- data.frame(bam = list.files(bam_dir, pattern = '.bam$', full.names = T)) %>%
  mutate(sample = sub("\\..*", "",basename(bam)))

peak_files <- data.frame(Peaks = list.files(peaks_dir, pattern = 'Peak$', full.names = T)) %>%
  mutate(SampleID = sub("\\..*", "",basename(Peaks))) %>%
  mutate(SampleID = gsub('_peaks', '', SampleID))

coldata_for_diffbind <- coldata %>% 
  filter(!is.na(control) & control != '') %>%
  # select(-description) %>%
  dplyr::rename(ControlID = control, SampleID = sample) %>% 
  separate(SampleID, into = c('sample', 'Replicate'), remove = F, sep = '_REP') %>%
  mutate(peakCaller = peak_caller) 

samplesheet <- coldata_for_diffbind %>%
  left_join(bam_files %>% select(SampleID = sample, bamReads = bam), by = 'SampleID') %>%
  left_join(bam_files %>% select(ControlID = sample, bamControl = bam), by = 'ControlID') %>%
  left_join(peak_files, by = 'SampleID')

diffbind_obj <- dba(sampleSheet = samplesheet, scoreCol = 5)
diffbind_count <- dba.count(diffbind_obj, bUseSummarizeOverlaps = TRUE)
qsave(diffbind_count, 'diffbind_count.qs')

