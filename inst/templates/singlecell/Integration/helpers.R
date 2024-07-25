# some code showing integrating two seurat objects
library(Seurat)
library(tidyverse)
library(harmony)

# We assume exp1 and exp2 has a Group column with naming the sample groups
# We create a batch annotation for each batch
exp1=readRDS("data/exp1.rds")
exp1$batch="n10"
exp2=readRDS("data/exp2.rds")
exp2$batch="n6"

# Normalize ----
exp = SCTransform(exp, verbose = FALSE,conserve.memory=TRUE,
                  variable.features.n = 3000)
exp <- RunPCA(exp)
ElbowPlot(exp,ndims=40)
end_dimension=35
resolution=0.5

# Clustering for each batch ----
exp <- FindNeighbors(exp, dims = 1:end_dimension, reduction = "pca")
exp <- FindClusters(exp, resolution = resolution, cluster.name = "unintegrated_clusters")
exp <- RunUMAP(exp, dims = 1:end_dimension, reduction = "pca", reduction.name = "umap.unintegrated")
saveRDS(exp, file="data/merged_umap.rds")

## Plot by batch----
DimPlot(exp, reduction = "umap.unintegrated", split.by = c("Group"),
        group.by = c("batch"))

# Integration ----
exp <- IntegrateLayers(
  object = exp, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
exp[["RNA"]] <- JoinLayers(exp[["RNA"]])
end_dimension=35
resolution=0.5
exp <- FindNeighbors(exp, reduction = "harmony", dims = 1:end_dimension)
exp <- FindClusters(exp, resolution = resolution, cluster.name = "harmony_clusters")
exp <- RunUMAP(exp, reduction = "harmony", dims = 1:end_dimension, reduction.name = "umap.harmony")
saveRDS(exp, file="data/integrated_harmony.rds")

## Plot by Group and Cluster ----
DimPlot(
  exp,
  reduction = "umap.harmony",
  split.by = c("Group"),
  group.by = c("batch"),
  combine = TRUE, label.size = 2
)

DimPlot(
  exp,
  reduction = "umap.harmony",
  split.by = c("Group"),
  group.by = c("ident"),
  combine = TRUE, label.size = 2
)

