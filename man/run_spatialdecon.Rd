% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolution_methods.r
\name{run_spatialdecon}
\alias{run_spatialdecon}
\title{Run SpatialDecon deconvolution method}
\usage{
run_spatialdecon(scrna_path, spatial_path, celltype_final, output_path)
}
\arguments{
\item{scrna_path}{Path to scRNA-seq data}

\item{spatial_path}{Path to spatial transcriptomics data}

\item{celltype_final}{Column name in sc_obj@meta.data containing cell type information}

\item{output_path}{Path to save output}
}
\value{
SpatialDecon result
}
\description{
Run SpatialDecon deconvolution method
}

\examples{
scrna_path = './sc_rna.h5seurat'
spatial_path = './st_data.h5seurat'
celltype_final = 'celltype'
output_path = 'output/'

sc_obj <- LoadH5Seurat(scrna_path)
spatial_obj <- LoadH5Seurat(spatial_path)
sc_count <- as.matrix(sc_obj@assays$RNA@counts)
spatial_count <- as.matrix(spatial_obj@assays$RNA@counts)

cellType <- sc_obj@meta.data[, celltype_final, drop = FALSE]
st_counts_norm = sweep(spatial_count, 2, colSums(spatial_count), "/") * mean(colSums(spatial_count))
st_object=CreateSeuratObject(counts=st_counts_norm,assay="Spatial")
stopifnot(setequal(colnames(st_object),rownames(spatial_location)))
st_object=AddMetaData(st_object,spatial_location[colnames(st_object),1],col.name="x")
st_object=AddMetaData(st_object,spatial_location[colnames(st_object),2],col.name="y")

stopifnot(all(colnames(sc_count)==rownames(cellType)))

sc_counts_matrix=as.matrix(sc_count)
sc_counts_matrix=Matrix::Matrix((sc_counts_matrix),sparse=TRUE)
sc_labels_df=data.frame(cell_barcodes=rownames(cellType),sc_labels=cellType$celltype)
sc_matrix <- create_profile_matrix(
    mtx = sc_counts_matrix,          
    cellAnnots = sc_labels_df,  
    cellTypeCol = "sc_labels",  
    cellNameCol = "cell_barcodes",           
    matrixName = "custom_cell_type_matrix", 
    outDir = NULL, 
    normalize = TRUE,
    minCellNum = 1,
    minGenes = 1
) 
  
    
## deconvolution
res = runspatialdecon(object = st_object,
                  bg = 0.01,
                  X = sc_matrix,
                  align_genes = TRUE)
weights=t(res$beta)
norm_weights=sweep(weights, 1, rowSums(weights), "/")
SpatialDecon_result <- norm_weights
save(SpatialDecon_result, file = file.path(output_path, "SpatialDecon_result.Rdata"))
}
