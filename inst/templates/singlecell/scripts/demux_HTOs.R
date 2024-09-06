library(Seurat)
library(qs)

############################ perform demultiplexing ############################

# replace these paths with ones to your data (cellranger outputs)
sample_matrix <- ReadMtx('data/tumor_4068_ovarian_cdki/processed_data/matrix.mtx.gz',
                         cells = 'data/tumor_4068_ovarian_cdki/processed_data/barcodes.tsv.gz',
                         features = 'data/tumor_4068_ovarian_cdki/processed_data/features.tsv.gz')

# create two matrices of counts: one of for hashtag oligo counts, and one for counts for actual genes
# include in the HTO count matrix only those HTOs that are actually assigned to samples in your dataset (in this case, Hashtag1 and Hashtag2)
hto_matrix <- full_matrix[grepl('Hashtag[12]+', rownames(full_matrix)), ]
expression_matrix <- full_matrix[!grepl('Hashtag', rownames(full_matrix)), ]

# create a Seurat object from the raw data, including a slot for HTO counts
sample_seurat <- CreateSeuratObject(
  counts = Matrix::Matrix(as.matrix(expression_matrix), sparse = T))
sample_seurat[["HTO"]] <- CreateAssayObject(counts = hto_matrix)

# normalize both slots in the Seurat object
sample_seurat <- NormalizeData(sample_seurat)
sample_seurat <- NormalizeData(sample_seurat, assay = "HTO", normalization.method = "CLR")

# perform demultiplexing. adjust positive.quantile as necessary to call more/fewer cells as hashtag-positive
sample_seurat <- HTODemux(sample_seurat, assay = "HTO", positive.quantile = 0.99)

qsave(sample_seurat, 'data/processed/hto_demux_seurat.qs')
# saveRDS(sample_seurat, 'data/processed/hto_demux_seurat.rds')


################## evaluate demultiplexing performance #########################

# distributions of expression of hashtags should make sense considering hashtag assigned
RidgePlot(sample_seurat, assay = "HTO", features = c("Hashtag1", 'Hashtag2'), ncol = 2) 

# evaluate expression of hashtags vs calls for singlet, doublet, and unassigned
FeatureScatter(sample_seurat, feature1 = "Hashtag1", feature2 = "Hashtag2")
HTOHeatmap(sample_seurat, assay = "HTO")

# evaluate nCount_RNA of cells classified as doublets
Idents(sample_seurat) <- "HTO_classification.global"
VlnPlot(sample_seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

