library(rstudioapi)
library(tidyverse)
setwd(fs::path_dir(getSourceEditorContext()$path))
# NOTE: This code will check version, this is our recommendation, it may work
#.      other versions
stopifnot(R.version$major>= 4) # requires R4
stopifnot(compareVersion(R.version$minor,"3.3")>=0) # requires >=4.3.3
stopifnot(compareVersion(as.character(BiocManager::version()), "3.18")>=0)
stopifnot(compareVersion(as.character(packageVersion("Seurat")), "5.0.0")>=0)

library(Seurat)
library(data.table)
library(hdf5r)

### Set up run information
data_dir <- "/path/to/cellranger/output/folders/"

samples <- c("sample1",  "sample2",  "sample3")

### Make individual seurat objects for each sample

for (i in 1:length(samples)){
  seurat_data <- Read10X_h5(paste(c(data_dir,samples[i],"/outs/raw_feature_bc_matrix.h5"),sep="",collapse = ""))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,  ## only keep cells with at least 100 genes
                                   project = samples[i])
  assign(paste0(samples[i], "_seurat"),
         seurat_obj) # stores Seurat object in variable of corresponding sample name
}

### Merge all seurat objects

seurat_ID <- paste0(samples, "_seurat") # get names of all objects


u <- get(seurat_ID[2])
for (i in 3:length(seurat_ID)) {
  u <- c(u, get(seurat_ID[i]))
} ## makes a list of all seurat objects

seurat_merge <- merge(x = get(seurat_ID[1]),
                      y = u,
                      add.cell.ids = all_samples,
                      project = "my_scRNA_project")


# Mitochondrial genes for mouse genome
idx <- grep("^mt-", rownames(GetAssay(seurat_merge, "RNA")))
rownames(GetAssay(seurat_merge, "RNA"))[idx]
# Mitochondrial genes vs. nuclear genes ratio
seurat_merge$mitoRatio <- PercentageFeatureSet(object = seurat_merge, pattern = "^mt-")
seurat_merge$mitoRatio <- seurat_merge@meta.data$mitoRatio/100 # Divide by 100 for Ratio instead of Percentage

# Number of genes per UMI for each cell
seurat_merge$Log10GenesPerUMI <- log10(seurat_merge$nFeature_RNA) / log10(seurat_merge$nCount_RNA)

# Extract cell level metadata
metadata <- seurat_merge@meta.data
metadata$barcode <- rownames(metadata)

# Import experimental metadata
metaexp <- read.csv("/path/to/experimental/metadata/meta.csv")

# Check matching of IDs
all(metaexp$sample %in% metadata$orig.ident)
all(metadata$orig.ident %in% metaexp$sample)

#change headings to match
colnames(metaexp)[1] <- "orig.ident"

metafull <- plyr::join(metadata, metaexp,
                       by = c("orig.ident"))

# Replace seurat object metadata
if(all(metafull$barcode == rownames(seurat_merge@meta.data))) {
  rownames(metafull) <- metafull$barcode
  seurat_merge@meta.data <- metafull
}


## Join layers (each sample is a separate layer)
seurat_merge[["RNA"]] <- JoinLayers(seurat_merge[["RNA"]])

### Save Seurat object for future processing
saveRDS(seurat_merge, file = "seurat_pre-filtered.rds")
write.csv(seurat_merge@meta.data, file = "metadata_pre-filtered.csv")
