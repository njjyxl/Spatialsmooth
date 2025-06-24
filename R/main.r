#' Run all deconvolution methods
#'
#' @import Matrix Seurat SeuratDisk reticulate
#' @param scrna_path Path to scRNA-seq data
#' @param spatial_path Path to spatial transcriptomics data
#' @param python_path Path to Python executable
#' @param output_path Directory to save output, If not set, save to the temporary directory tempdir().
#' @param celltype_final Column name in sc_obj@meta.data containing cell type information
#' @return all deconvolution result
#' @export
run_deconvolution_methods <- function(scrna_path, spatial_path, celltype_final, python_path, output_path = tempdir()) {
  if (!dir.exists(output_path)) {
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Define method configurations
  methods <- list(
    CARD = list(func = run_card, use_python = FALSE),
    RCTD = list(func = run_rctd, use_python = FALSE),
    SPOTlight = list(func = run_spotlight, use_python = FALSE),
    SpatialDWLS = list(func = run_spatialdwls, use_python = TRUE),
    AdRoit = list(func = run_adroit, use_python = FALSE),
    Seurat = list(func = run_seurat, use_python = FALSE),
    Redeconve = list(func = run_redeconve, use_python = FALSE),
    SpatialDecon = list(func = run_spatialdecon, use_python = FALSE)
  )
  
  # A new way to get registered
  if (exists("registered_methods", envir = .GlobalEnv)) {
    reg_methods <- get("registered_methods", envir = .GlobalEnv)
    for (method_name in names(reg_methods)) {
      m_info <- reg_methods[[method_name]]
      if (exists(m_info$function_name, envir = .GlobalEnv)) {
        methods[[method_name]] <- list(
          func = get(m_info$function_name, envir = .GlobalEnv),
          use_python = m_info$use_python
        )
      }
    }
  }
  
  # Worker function to run a single method
  run_worker <- function(method_name, method_info) {
    tryCatch({
      cat(paste("Running", method_name, "...\n"))
      if (method_info$use_python) {
        method_info$func(scrna_path, spatial_path, celltype_final, python_path, output_path)
      } else {
        method_info$func(scrna_path, spatial_path, celltype_final, output_path)
      }
      cat(paste(method_name, "completed successfully.\n"))
      return(list(name = method_name, status = "success"))
    }, error = function(e) {
      msg <- paste("Error in", method_name, ":", conditionMessage(e))
      cat(msg, "\n")
      return(list(name = method_name, status = "error", error = msg))
    })
  }
  
  # Setup parallel execution
  suppressMessages(library(doParallel))
  suppressMessages(library(foreach))
  suppressMessages(library(parallel))
  num_cores <- min(length(methods), detectCores() - 1)  # Leave one core free
  cl <- makeCluster(num_cores)
    
  # Get function name
  all_func_names <- unique(c(
    "run_card", "run_rctd", "run_spotlight", "run_spatialdwls",
    "run_adroit", "run_seurat", "run_redeconve", "run_spatialdecon"
  ))
  
  # Add the function name of the registration method
  if (exists("registered_methods", envir = .GlobalEnv)) {
    reg_methods <- get("registered_methods", envir = .GlobalEnv)
    custom_funcs <- sapply(reg_methods, function(m) m$function_name)
    all_func_names <- unique(c(all_func_names, custom_funcs))
  }
    
  clusterEvalQ(cl, {
    suppressMessages(library(data.table))
    suppressMessages(library(Seurat))
    suppressMessages(library(SeuratDisk))
    suppressMessages(library(CARD))
    suppressMessages(library(SPOTlight))
    suppressMessages(library(Giotto))
    suppressMessages(library(Redeconve))
    suppressMessages(library(AdRoit))
    suppressMessages(library(SpatialDecon))
    # Load the additional package specified by the registration method.
    if (exists("registered_methods", envir = .GlobalEnv)) {
      reg_methods <- get("registered_methods", envir = .GlobalEnv)
      for (pkg in unique(unlist(lapply(reg_methods, function(m) m$packages)))) {
        suppressMessages(library(pkg, character.only = TRUE))
      }
    }
    NULL
  })
                                
  # Export all functions
  clusterExport(cl, varlist = c(
    "scrna_path", "spatial_path", "celltype_final", 
    "python_path", "output_path", all_func_names
  ), envir = environment())
  
  # Execute methods in parallel
  results <- parLapplyLB(cl, names(methods), function(method_name) {
    run_worker(method_name, methods[[method_name]])
  })
  
  stopCluster(cl)
  
  # Print summary
  cat("\n===== Execution Summary =====\n")
  for (res in results) {
    if (res$status == "success") {
      cat(paste("[SUCCESS]", res$name, "\n"))
    } else {
      cat(paste("[FAILED]", res$name, "->", res$error, "\n"))
    }
  }
  cat("=============================\n")
}
