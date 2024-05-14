# prepare seurat object with all 3 modes

# read cellranger-arc outs
#arc_h5 <- '../aspera_download_sep23/SN0288408/KW11249_atsibris/230628_10X_KW11252_bcl/cellranger-arc-2.0.0/GRCh38/BRI-2286/outs/filtered_feature_bc_matrix.h5'
arc_h5 <- '../final/batch3_atac/cellranger-arc_output/outs/filtered_feature_bc_matrix.h5'

arc_data <- Read10X_h5(arc_h5)

# read cellranger outs
#h5path <- '../aspera_download_sep23/SN0288408/KW11249_atsibris/230711_10X_KW11249_bcl/cellranger-7.1.0/GRCh38/BRI-2285_hashing/outs/filtered_feature_bc_matrix.h5'
h5path <- '../final/batch3_adt/AtheTeaSeqMulti/outs/multi/count/raw_feature_bc_matrix.h5'
counts <- Read10X_h5(h5path)

# select cell barcodes detected by both ATAC and ADT
joint.bcs <- intersect(colnames(arc_data$Peaks), colnames(counts$`Antibody Capture`))

so <- CreateSeuratObject(counts = arc_data$`Gene Expression`[, joint.bcs])

# add cellranger ATAC
# fragpath <- '../aspera_download_sep23/SN0288408/KW11249_atsibris/230628_10X_KW11252_bcl/cellranger-arc-2.0.0/GRCh38/BRI-2286/outs/atac_fragments.tsv.gz'
fragpath <- '../final/batch3_atac/cellranger-arc_output/outs/atac_fragments.tsv.gz'
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

so[["CellRangerPeaks"]] <- CreateChromatinAssay(
  counts = arc_data$Peaks[, joint.bcs],
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# add HTO
so[["HTO"]] <- CreateAssayObject(counts = counts$`Multiplexing Capture`[c('Hashtag1', 'Hashtag2','Hashtag3', 'Hashtag4', 'Hashtag5','Hashtag6','Hashtag7'), joint.bcs])
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
so <- NormalizeData(so, assay = "HTO", normalization.method = "CLR")
so <- HTODemux(so, assay = "HTO", positive.quantile = 0.99)

# table(so@meta.data$orig.ident, so@meta.data$hash.ID) %>% as.matrix() %>% kable() %>% kable_styling()

# add ADT data
so[["ADT"]] <- CreateAssayObject(counts = counts$`Antibody Capture`[, joint.bcs])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR", margin = 2)

qsave(so, "data/tea_seurat.qs")

hashtag_counts <- table(so@meta.data$hash.ID) %>% as.data.frame() %>%
  set_colnames(c('Hashtag', 'Number of cells')) %>% dplyr::filter(grepl("Hash", Hashtag))
write.csv(hashtag_counts, 'data/hashtag_counts.csv', row.names = FALSE)

