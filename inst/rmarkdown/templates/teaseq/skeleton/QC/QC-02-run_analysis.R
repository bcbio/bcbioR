library(Seurat)
library(tidyverse)
library(Matrix)
library(data.table)
library(magrittr)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(qs)
library(bcbioR)
# Run analysis on multi-omic sc Data
seurat <- qread("data/tea_seurat.qs")

# keep singlet only
seurat <- subset(x = seurat, subset = HTO_classification.global == "Singlet")
# table(seurat@meta.data$hash.ID)

# add sampleinfo to metadata
seurat@meta.data <- seurat@meta.data %>%
  rownames_to_column("cellID") %>%
  left_join(sampleinfo, by = "hash.ID") %>%
  column_to_rownames("cellID")

#Compute percent mito ratio
seurat$percent_mito <- PercentageFeatureSet(seurat, pattern = "^MT")
# Compute total reads
seurat$total_reads <- Matrix::colSums(GetAssayData(seurat, slot = "counts"))
# Compute novelty score
seurat$novelty <- log10(seurat@meta.data$nFeature_RNA)/log10(seurat@meta.data$nCount_RNA)

seurat <- SCTransform(seurat, conserve.memory=TRUE)

# cell cycle
cc_genes <- cc.genes.updated.2019

# Extract phase specific genes
s.genes <- cc_genes$s.genes
g2m.genes <- cc_genes$g2m.genes

# Perform cell cycle scoring
seurat <- CellCycleScoring(object = seurat,
                           s.features = s.genes,
                           g2m.features = g2m.genes,
                           set.ident = TRUE)

# clustering
seurat <- RunPCA(object = seurat, verbose = T)
pct <- seurat@reductions$pca@stdev/sum(seurat@reductions$pca@stdev) * 100
co1 <- which(cumsum(pct) > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# co <- min(co1, co2)
co <- 30
seurat <- RunPCA(object = seurat, npcs = co, verbose = TRUE)
seurat <- FindNeighbors(object = seurat, reduction = "pca", dims = 1:co)
seurat <- RunUMAP(object = seurat, dims=1:co )

# clustering by ADT data # Error: didn't do CLR norm because copied the code from 1st dataset
seurat <- ScaleData(seurat, assay = "ADT")
seurat <- RunPCA(seurat, features = rownames(seurat@assays$ADT@counts), verbose = TRUE, assay= "ADT", reduction.name = "pcaADT", reduction.key = "pcaADT_")
seurat <- RunUMAP(seurat, dims = 1:20, reduction = "pcaADT", assay = "ADT", reduction.name = "umapADT", reduction.key = "umapADT_")

# calling peaks using MACS2, following James's code
DefaultAssay(seurat)<-"CellRangerPeaks"
peaks <- CallPeaks(seurat) # took 16 min
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
# peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)


# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seurat),
  features = peaks,
  cells = colnames(seurat)
)

#qsave(seurat,"data/tmp.qs")
# ATAC data normalization, dimension reduction, and clustering
DefaultAssay(seurat)<-"MACS2peaks"
seurat <- RunTFIDF(seurat, assay = "MACS2peaks")
seurat <- FindTopFeatures(seurat, min.cutoff = 'q0', assay = "MACS2peaks")
seurat <- RunSVD(seurat, assay = "MACS2peaks")
seurat <- RunUMAP(object = seurat, dims = 2:30, reduction = 'lsi', assay = "MACS2peaks", reduction.name = "umapATAC", reduction.key = "umapATAC_" )
seurat <- FindNeighbors(object = seurat, reduction = 'lsi', dims = 2:30)
seurat <- FindClusters(object = seurat, verbose = FALSE, algorithm = 3)

seurat$blacklist_fraction <- FractionCountsInRegion(
  object = seurat,
  assay = 'MACS2peaks',
  regions = blacklist_hg38_unified
)

## add ATAC per barcode metrics
pbm <- read.csv(
  file = "cellranger-arc_output/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

seurat@meta.data <- seurat@meta.data %>%
  rownames_to_column("cellID") %>%
  left_join(pbm %>% rownames_to_column('cellID'), by = "cellID") %>%
  column_to_rownames("cellID")

seurat$pct_reads_in_peaks <- seurat$atac_peak_region_fragments / seurat$atac_fragments * 100

# compute nucleosome signal score per cell
seurat <- NucleosomeSignal(object = seurat) # added "nucleosome_signal", "nucleosome_percentile" to metadata
# compute TSS enrichment score per cell
seurat <- TSSEnrichment(object = seurat, fast = FALSE) # added "TSS.enrichment", "TSS.percentile" to metadata
qsave(seurat, 'data/tea_seurat_unfiltered_clustered.qs')
