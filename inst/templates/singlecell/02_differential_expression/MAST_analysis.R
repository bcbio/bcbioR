##### library loading #####
library(Seurat)
library(apeglm)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(MAST)
library(tidyverse)
library(data.table)
library(lme4)
library(R.utils)
library(glue)
library(optparse)
##### parameter parse #####
options(stringsAsFactors = F)
option_list = list(make_option("--seurat_obj", default = "https://github.com/bcbio/bcbioR-test-data/raw/refs/heads/main/singlecell/tiny.rds"),
                   make_option("--resolution_column", default = "integrated_snn_res.0.4"),
                   make_option("--cluster_name", default = "2"),
                   make_option("--contrast", default = "age"),
                   make_option("--outputDir", default = "out")
)
args = parse_args(OptionParser(option_list = option_list))

invisible(list2env(args,environment()))
column <- contrast
system(glue("mkdir -p {outputDir}"))


message("[Preparing inputs for MAST modeling]")
##### Read in provided seurat #####
if (isUrl(seurat_obj)){
  seurat <- readRDS(url(seurat_obj))  
}else{
  seurat <- readRDS(seurat_obj)
}

message("Input seurat object: ",seurat_obj)
DefaultAssay(seurat) <- "RNA"
message("RNA is set as the default assay")
message("Column name of clustering to use: ",resolution_column)
Idents(object = seurat) <- resolution_column
message("Subset original seurat to be only cluster ",cluster_name," for faster computing!")
data_subset <- subset(x = seurat, idents = cluster_name)
existintLayers <- Layers(data_subset[["RNA"]])
if(existintLayers!="counts"){
  print("Not only counts excisted as layers in this object")
  print("Make sure your default slot is counts and it is raw counts")
}
##### Start from raw count for MAST #####
message("Natural log of raw counts with pseudobulk 1 used for MAST modeling")
sce <- as.SingleCellExperiment(data_subset) 
##### Log-Normalize Seurat for visualization later #####
message("Total counts normalization and log1p transformation done to raw counts")
message("New layer Data added to the seurat object for visualization later!")
data_subset <- NormalizeData(
  data_subset,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1)
saveRDS(data_subset,glue("{outputDir}/processed_seurat.rds"))
##### Continue MAST input prep #####
assay(sce, "log") = log(counts(sce) + 1)
# Scaling ngenes
cdr = colSums(assay(sce, "log")>0)
colData(sce)$cngeneson = scale(cdr)

# Create new sce object (only 'log' count data)
sce.1 = SingleCellExperiment(assays = list(log = assay(sce, "log")))
colData(sce.1) = colData(sce)

#change to sca
sca = SceToSingleCellAssay(sce.1)

message("Subset genes observed in at least 10% of cells")

expressed_genes <- freq(sca) > 0.1
sca_filtered <- sca[expressed_genes, ] 

cdr2 <- colSums(SummarizedExperiment::assay(sca_filtered)>0)

SummarizedExperiment::colData(sca_filtered)$ngeneson <- scale(cdr2)
SummarizedExperiment::colData(sca_filtered)$orig.ident <-
  factor(SummarizedExperiment::colData(sca_filtered)$orig.ident)
SummarizedExperiment::colData(sca_filtered)[[column]] <-
  factor(SummarizedExperiment::colData(sca_filtered)[[column]])

##### MAST modeling #####
message("[MAST modeling with supplied contrasts]")

message("Note: this step is time-consuming!")

comp_name <- levels(SummarizedExperiment::colData(sca_filtered)[[column]])[2]
lrt_name <- paste0(column, comp_name)
formula_touse <- as.formula(paste0("~ ngeneson + (1 | orig.ident) + ", column))

zlmCond <- suppressMessages(MAST::zlm(formula_touse,  sca_filtered, method='glmer',
                                      ebayes = F,strictConvergence = FALSE))
summaryCond_column <- suppressMessages(MAST::summary(zlmCond,doLRT=lrt_name))

##### MAST outputs #####
message("[Main MAST computation done, result outputs]")

summary_cond_file = paste0(outputDir,"/MAST_RESULTS_",cluster_name,"_", column, ".rds")
saveRDS(summaryCond_column, file = summary_cond_file)

message("Full MAST object saved to file ", summary_cond_file)


summaryDt_column <- summaryCond_column$datatable
fcHurdle_column <- merge(summaryDt_column[contrast == lrt_name & component == 'H',
                                          .(primerid, `Pr(>Chisq)`)], 
                         # This extracts hurdle p-values 
                         summaryDt_column[contrast == lrt_name & component == 'logFC', 
                                          .(primerid, coef, ci.hi, ci.lo)], 
                         # This extract LogFC data
                         by = 'primerid')
fcHurdle_column <- stats::na.omit(as.data.frame(fcHurdle_column))
fcHurdle_column$fdr <- p.adjust(fcHurdle_column$`Pr(>Chisq)`, 'fdr')
to_save_column <- fcHurdle_column

full_res_file = paste0(outputDir,"/FULL_MAST_RESULTS_",cluster_name,"_", column, ".csv")
write.table(to_save_column, file=full_res_file, row.names = FALSE, sep=",")


message("MAST summary results output to csv files")


fcHurdleSig_column <- merge(fcHurdle_column[fcHurdle_column$fdr < .05,],
                            as.data.table(mcols(sca_filtered)),
                            by = 'primerid')

setorder(fcHurdleSig_column, fdr)

sig_res_file = paste0(outputDir,"/SIG_MAST_RESULTS_padj<0.05_",cluster_name,"_", column, ".csv")
write.table(fcHurdleSig_column, file=sig_res_file, row.names = FALSE, sep=",")

message("Significant MAST summary results output to csv files")
