#' Run all deconvolution methods
#'
#' @import Matrix Seurat SeuratDisk reticulate
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param python_path Path to Python executable
#' @param output_path Directory to save output
#' @param celltype_final Column name in sc_obj@meta.data containing cell type information
#' @return all deconvolution result
#' @export
run_deconvolution_methods <- function(scrna_path, spatial_path, celltype_final, python_path, output_path) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  methods <- list(
    CARD = run_card,
    RCTD = run_rctd,
    SPOTlight = run_spotlight,
    SpatialDWLS = run_spatialdwls,
    AdRoit = run_adroit,
    Seurat = run_seurat,
    Redeconve = run_redeconve,
    SpatialDecon = run_spatialdecon
  )
  
  for (method_name in names(methods)) {
    tryCatch({
      cat(paste("Running", method_name, "...\n"))
      method_func <- methods[[method_name]]
      
      if (method_name == "SpatialDWLS") {
        method_func(scrna_path, spatial_path, celltype_final, python_path, output_path)
      } else {
        method_func(scrna_path, spatial_path, celltype_final, output_path)
      }
      
      cat(paste(method_name, "completed successfully.\n"))
    }, error = function(e) {
      cat(paste("Error in running", method_name, ":", conditionMessage(e), "\n"))
    })
  }
  
  cat("All methods completed.\n")
}

