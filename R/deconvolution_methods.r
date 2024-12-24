#' Run CARD deconvolution method
#'
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param output_path Path to save output
#'
#' @return CARD object
#' @export
run_card <- function(scrna_path, spatial_path, celltype_final, output_path) {
  sc_obj <- LoadH5Seurat(scrna_path)
  spatial_obj <- LoadH5Seurat(spatial_path)
  sc_count <- as.matrix(sc_obj@assays$RNA@counts)
  spatial_count <- as.matrix(spatial_obj@assays$RNA@counts)
  spatial_location <- data.frame(spatial_obj@meta.data)
  colnames(spatial_location) <- c('x', 'y')
  cellType <- sc_obj@meta.data

  CARD_obj = createCARDObject(
    sc_count = sc_count,
    sc_meta = cellType,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "celltype",
    ct.select = unique(cellType$celltype),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5
  )
  
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  CARD_result <- CARD_obj@Proportion_CARD
  save(CARD_result, file = file.path(output_path, "CARD_result.Rdata"))
  return(CARD_obj)
}

#' Run RCTD deconvolution method
#'
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param output_path Path to save output
#'
#' @return RCTD result
#' @export
run_rctd <- function(scrna_path, spatial_path, celltype_final, output_path) {
  sc_obj <- LoadH5Seurat(scrna_path)
  diff_list <- list()
  c = 0
  for (i in seq_along(unique(sc_obj@meta.data[,celltype_final]))){
      if(sum(sc_obj@meta.data[,celltype_final] == unique(sc_obj@meta.data[,celltype_final])[i]) < 1){
          c = c+1
          diff_list[[c]] <- unique(sc_obj@meta.data[,celltype_final])[i]
          print(unique(sc_obj@meta.data[,celltype_final])[i])
      }
  }
  for (i in diff_list){
          sc_obj = sc_obj[,sc_obj@meta.data[,celltype_final]!=i]
      }
  sc_obj@meta.data[,celltype_final] <- as.factor(as.character(sc_obj@meta.data[,celltype_final]))
  ### Load in/preprocess your data, this might vary based on your file type
  print('prepare data')
  counts <- data.frame(sc_obj@assays$RNA@counts)
  colnames(counts) <- colnames(sc_obj)
  meta_data <- data.frame(sc_obj@meta.data[,celltype_final])
  cell_types <- meta_data[,1]
  # 使用 gsub 函数替换或移除特殊字符
  cell_types <- as.character(cell_types)
  cell_types <- gsub("/", "_", cell_types) # 将 / 替换为 _
  # cell_types <- as.factor(cell_types)
  names(cell_types) <- rownames(sc_obj@meta.data)
  cell_types <- as.factor(cell_types)
#   nUMI_df <- data.frame(colSums(sc_obj@assays$RNA@counts))
  nUMI_df <- data.frame(colSums(counts))
#   nUMI <- nUMI_df$colSums.sc_obj.assays.RNA
  nUMI <- nUMI_df$colSums
  names(nUMI) <- rownames(nUMI_df)

  ### Create the Reference object
  # reference <- Reference(counts, cell_types, nUMI, min_UMI=1)
  reference <- Reference(counts, cell_types, nUMI)
  #> Warning in Reference(counts, cell_types, nUMI): Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
  #>             is intended, there is no problem.
  spatial_obj <- LoadH5Seurat(spatial_path)
  counts <- data.frame(spatial_obj@assays$RNA@counts) # load in counts matrix
  colnames(counts) <- colnames(spatial_obj)
  coords <- data.frame(colnames(spatial_obj))
  coords <- data.frame(spatial_obj@meta.data)
  nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
  
  puck <- SpatialRNA(coords, counts, nUMI)
  myRCTD <- create.RCTD(puck, reference, max_cores = 8, CELL_MIN_INSTANCE = 1)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  weights <- as.matrix(results$weights)
  # normalize the cell type proportions to sum to 1.
#   norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/')
  norm_weights = sweep(results$weights, 1, rowSums(weights), '/') 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  RCTD_result <- norm_weights
  
#   return(list(myRCTD = myRCTD, 
#               results = results, 
#               norm_weights = norm_weights, 
#               cell_type_names = cell_type_names, 
#               spatialRNA = spatialRNA))
  # 假设结果保存在RCTD_result变量中
  save(RCTD_result, file = file.path(output_path, "RCTD_result.Rdata"))
}


#' Run SPOTlight deconvolution method
#'
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param output_path Path to save output
#'
#' @return SPOTlight result
#' @export
run_spotlight <- function(scrna_path, spatial_path, celltype_final, output_path) {
  sc <- LoadH5Seurat(scrna_path)
  st <- LoadH5Seurat(spatial_path)

  set.seed(123)
  sc <- Seurat::SCTransform(sc, verbose = FALSE)

  Idents(sc) <- sc@meta.data[,celltype_final]

  cluster_markers_all <- FindAllMarkers(object = sc, 
                                                assay = "SCT",
                                                slot = "data",
                                                verbose = TRUE, 
                                                only.pos = TRUE)
  spatial_location <- data.frame(colnames(st))
  spatial_location <- data.frame(st@meta.data)
  spotlight_ls <- spotlight_deconvolution(
    se_sc = sc,
    counts_spatial = st@assays$RNA@counts,
    clust_vr = celltype_final,
    cluster_markers = cluster_markers_all,
    cl_n = 100,
    hvg = 3000,
    ntop = NULL,
    transf = "uv",
    method = "nsNMF",
    min_cont = 0
  )
  
  decon_mtrx <- spotlight_ls[[2]]
  rownames(decon_mtrx) <- colnames(st@assays$RNA@counts)
  SPOTlight_result <- decon_mtrx[, which(colnames(decon_mtrx) != "res_ss")]
#   return(list(spotlight_ls = spotlight_ls, SPOTlight_result = SPOTlight_result))
  # 假设结果保存在SPOTlight_result变量中
  save(SPOTlight_result, file = file.path(output_path, "SPOTlight_result.Rdata"))
}

#' Run SpatialDWLS deconvolution method
#'
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param python_path Path to Python executable
#' @param output_path Path to save output
#'
#' @return SpatialDWLS result
#' @export
run_spatialdwls <- function(scrna_path, spatial_path, celltype_final, python_path, output_path) {
  instrs = createGiottoInstructions(python_path = python_path)
  
  sc_obj <- LoadH5Seurat(scrna_path)
  spatial_obj <- LoadH5Seurat(spatial_path)
  st_data <- createGiottoObject(
      raw_exprs = spatial_obj@assays$RNA@counts,
      instructions = instrs
  )
  st_data <- normalizeGiotto(gobject = st_data)
  st_data <- calculateHVG(gobject = st_data)
  gene_metadata = fDataDT(st_data)
#   featgenes = gene_metadata[hvg == 'yes']$gene_ID
  # featgenes = gene_metadata[hvg == 'no']$gene_ID
  gene_metadatas <- data.frame(gene_metadata)
  featgenes <- gene_metadatas[gene_metadatas$hvg == 'yes', "gene_ID"]

  st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
  st_data <- runUMAP(st_data, dimensions_to_use = 1:10)
  st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
  st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)

  sc_data <- createGiottoObject(
      raw_exprs = sc_obj@assays$RNA@counts,
      instructions = instrs
  )
  sc_data <- normalizeGiotto(gobject = sc_data)
  sc_data <- calculateHVG(gobject = sc_data)
  gene_metadata = fDataDT(sc_data)
#   featgenes = gene_metadata[hvg == 'yes']$gene_ID
  # featgenes = gene_metadata[hvg == 'no']$gene_ID
  gene_metadatas <- data.frame(gene_metadata)
  featgenes <- gene_metadatas[gene_metadatas$hvg == 'yes', "gene_ID"]
  sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
  sc_data@cell_metadata$leiden_clus <- as.character(sc_obj@meta.data[,celltype_final])
  
  gini_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
                                                     method = 'gini',
                                                     expression_values = 'normalized',
                                                     cluster_column = 'leiden_clus',
                                                     min_genes = 20,
                                                     min_expr_gini_score = 0.5,
                                                     min_det_gini_score = 0.5)
#   topgenes_gini = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']
  topgenes_gini <- gini_markers_subclusters %>%
  group_by(cluster) %>%
  slice_head(n = 100) %>%
  ungroup()
  topgenes_gini <- data.frame(topgenes_gini)
  sc_norm_exp <- 2^(sc_data@norm_expr)-1
  ExprSubset<-sc_norm_exp[as.character(topgenes_gini$genes),]
  Sig<-NULL
  for (i in as.character(unique(sc_obj@meta.data[,celltype_final]))){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(sc_obj@meta.data[,celltype_final]==i)]))))
  }
#   print(dim(ExprSubset))
#   print(dim(Sig))
#   print(length(unique(sc_obj@meta.data[,celltype_final])))

#   Sig <- as.matrix(Sig)
  colnames(Sig)<-as.character(unique(sc_obj@meta.data[,celltype_final]))
                          
  st_data <- runDWLSDeconv(st_data,sign_matrix = Sig, n_cell = 20)
  # 提取 DWLS 数据
  dwls_data <- as.matrix(st_data@spatial_enrichment$DWLS)
  # 将第一列设为行名
  rownames(dwls_data) <- dwls_data[, 1]
  # 删除第一列
  SpatialDWLS_result <- dwls_data[, -1]
  # 假设结果保存在SpatialDWLS_result变量中
  save(SpatialDWLS_result, file = file.path(output_path, "SpatialDWLS_result.Rdata"))
}
#' Run AdRoit deconvolution method
#'
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param output_path Path to save output
#'
#' @return AdRoit result
#' @export
run_adroit <- function(scrna_path, spatial_path, celltype_final, output_path) {
    sc_obj <- LoadH5Seurat(scrna_path)
    spatial_obj <- LoadH5Seurat(spatial_path)
    sc_count <- as.matrix(sc_obj@assays$RNA@counts)
    spatial_count <- as.matrix(spatial_obj@assays$RNA@counts)
    spatial_location <- data.frame(colnames(spatial_obj))
    spatial_location <- data.frame(spatial_obj@meta.data)
    
    colnames(spatial_location) <- c('x', 'y')
    cellType <-  sc_obj@meta.data[,celltype_final]
    intersect(rownames(spatial_count),rownames(sc_count))
    
    sc_count <- as.matrix(sc_count)
    spatial_count <- as.matrix(spatial_count)
    
    single.ref = ref.build(sc_count,
                       cellType,
                       genes=intersect(rownames(spatial_count),rownames(sc_count)),
                       normalize = "None",
                       multi.sample.bulk = T,
                       multi.sample.single = T,
                       nbootsids=5, 
                       minbootsize=50,
                       silent = F)

    print("single cell ref done")
    st_count_ma=as.matrix(spatial_count)
    st_count_ma = st_count_ma[intersect(rownames(st_count_ma),rownames(sc_count)),]
  
    registerDoSEQ()
    result=AdRoit.est(
         st_count_ma,
         single.ref,
         use.refvar = FALSE,
         per.sample.adapt = FALSE,
         silent = TRUE)

    results=t(result%>% as.data.frame(check.names=F))%>% as.data.frame(check.names=F)
    AdRoit_result <- results
    save(AdRoit_result,file = file.path(output_path,"AdRoit_result.Rdata"))
}

#' Run Seurat deconvolution method
#'
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param output_path Path to save output
#'
#' @return Seurat result
#' @export
run_seurat <- function(scrna_path, spatial_path, celltype_final, output_path) {
    sc_rna <- LoadH5Seurat(scrna_path)
    spatial <- LoadH5Seurat(spatial_path)
    
    sc_rna <- SCTransform(sc_rna)
    spatial <- SCTransform(spatial)
    anchors <- FindTransferAnchors(reference = sc_rna, 
                                   query = spatial, 
                                   k.anchor = 10, 
                                   k.filter = 200)
  
    predictions <- TransferData(anchorset = anchors, 
                                refdata = sc_rna@meta.data[,celltype_final], 
                                dims = 1:30)
    Seurat_result <- predictions[,which(colnames(predictions) != "predicted.id" & colnames(predictions) != "prediction.score.max")]
    colnames(Seurat_result) = gsub("prediction.score.", "", colnames(Seurat_result))
    save(Seurat_result, file = file.path(output_path,"Seurat_result.Rdata"))
}

#' Run Redeconve deconvolution method
#'
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param output_path Path to save output
#'
#' @return Redeconve result
#' @export
run_redeconve <- function(scrna_path, spatial_path, celltype_final, output_path) {
    sc_obj <- LoadH5Seurat(scrna_path)
    spatial_obj <- LoadH5Seurat(spatial_path)
    sc_count <- as.matrix(sc_obj@assays$RNA@counts)
    spatial_count <- as.matrix(spatial_obj@assays$RNA@counts)
    spatial_location <- data.frame(colnames(spatial_obj))
    spatial_location <- data.frame(spatial_obj@meta.data)
    
    colnames(spatial_location) <- c('x', 'y')
    cellType <-  sc_obj@meta.data
    cellType$barcodes <- rownames(cellType)
    cellType <- cellType[, c("barcodes", colnames(cellType)[-ncol(cellType)])]
    
    ## get reference
    ref <- get.ref(sc_count, cellType, dopar = F)
  
    # 执行CARD解卷积
    ## deconvolution
    res.ct <- deconvoluting(ref, spatial_count, genemode = "def", hpmode = "auto", dopar = T, ncores = 8)
    res.prop <- to.proportion(res.ct)
    Redeconve_result <- t(res.prop)
    save(Redeconve_result, file = file.path(output_path, "Redeconve_result.Rdata"))
}
