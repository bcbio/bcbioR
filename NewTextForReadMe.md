# SpatialViz

SpatialViz is a spatial single-cell RNA-seq visualization tool. In order to use the tool, you will need to clone this repository locally. When you launch the app it should look like:

<p align="center">
<img src="img/SpatialViz_default.png" width="600">
</p>

You can click on the <kbd>Browse...</kbd> button in order to upload your data. Two different sample data have been supplied in the "sample_data" directory. Each sample data RData object needs two data frames. Some example code for extracting these data frames from CosMx and Visium Seurat objects follows:

Example code extracting the two data frames from a CosMx Seurat object. 

xy_celltypes<-data.frame(matrix(ncol=4, nrow=length(colnames(SeuratObject)))#create an empty dataframe

                         
                  
metadata<-as.data.frame(SeuratObject@meta.data)#extract metadata from the seurat object

colnames(xy_celltypes)<-c("cellID", "x", "y", "celltype")


fill the dataframe

xy_celltypes[,1]<-colnames(SeuratObject)

xy_celltypes[,2]<-metadata$x_slide_coords

xy_celltypes[,3]<-metadata$y_slide_coords

xy_celltypes[,4]<-factor(metadata$CellAnnotation)#youre cell annotation column will have a different column name


remove any symbols from cell annotation names, for example:

levels(xy_celltypes$celltype) <- gsub("\\+", "", levels(xy_celltypes$celltype))#removing + sign from labels


expression_data <- as.data.frame(GetAssayData(SeuratObject, assay = "RNA", layer = "data"))#extract the expression data


save(xy_celltypes, expression_data, file="/path/to/save/sample_data.RData")#save the two dataframes as an RData object




Sample code to extract the two data frames from a Seurat Visium object


subset a single slide from the object, if there are multiple slides included in the object

if the slide names are in the orig.ident slot

Idents(SeuratObject)<-"orig.ident"


subSO<-subset(SeuratObject, idents="slideC")#your slide will have a different name


xy_celltypes<-data.frame(matrix(ncol=4, nrow=length(colnames(subSO))))#create an empty dataframe

the table should be called "xy_cellypes" and have four columns:
1. cellID
2. x
3. y
4. celltype

colnames(xy_celltypes)<-c("cellID", "x", "y", "celltype")



xy_celltypes[,1]<-colnames(subSO)

xy_celltypes[,2]<-subSO@images[["sliceC"]]@coordinates[["imagerow"]]#your images slot may have a different name

xy_celltypes[,3]<-subSO@images[["sliceC"]]@coordinates[["imagecol"]]#your images slot may have a different name


xy_celltypes[,4]<-subSO$predictionsID#your cell annotations factor column will have a different name

make sure there are no symbols or spaces in cell annotation names


extract expression data

expression_data <- as.data.frame(GetAssayData(subSO, assay = "Spatial", layer = "data"))#extract the expression data


save(expression_data, xy_celltypes, file="/path/to/save/sample_data.RData")#save the two data frames as an .RData object.


One table should be called "xy_cellypes" and have four columns:
1. cellID
2. x
3. y
4. celltype

> Note: It is important that the names of cell types lack an symbols like "+" or "-" as R can have trouble with these symbols.

The second data frame should be called "expression_data" and have rows for each gene and columns for each cell. The cellID from "xy_cellypes" should match the column names in "expression_data".



## Contributors

- Will Gammerdinger
- Lorena Pantano
- James Billingsley

## License

The code is freely available under the [MIT license](https://opensource.org/licenses/mit-license.html).
