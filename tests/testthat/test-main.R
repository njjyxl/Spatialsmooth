Sys.setenv(RETICULATE_PYTHON = "D:/Application/anaconda3/envs/python_3.10")
options(warn = -1)
suppressMessages(library(tidyverse))
suppressMessages(library(spacexr))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(MuSiC))
suppressMessages(library(CARD))
suppressMessages(library(SPOTlight))
suppressMessages(library(Giotto))
suppressMessages(library(Redeconve))
suppressMessages(library(AdRoit))
suppressMessages(library(SpatialDecon))
suppressMessages(library(gstat))
suppressMessages(library(dplyr))
suppressMessages(library(future))
suppressMessages(library(spatialreg))
suppressMessages(library(sp))
suppressMessages(library(igraph))
suppressMessages(library(fields))
suppressMessages(library(reticulate))

suppressMessages(library(Spatialsmooth))

test_that("deconvolution_module", {
  scrna_path <- './sc_rna.h5seurat'
  spatial_path <- './st_data.h5seurat'

  skip_if_not(file.exists(scrna_path), "scRNA data file not found")
  skip_if_not(file.exists(spatial_path), "Spatial data file not found")

  # 获取当前工作目录
  current_dir <- getwd()
  # 构建输出路径
  output_path <- file.path(current_dir, "test_output")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  # 使用系统上可用的Python路径
  # 在实际测试中，您可能需要设置一个可用的路径或使用条件跳过
  python_path <- Sys.which("D:/Application/anaconda3/envs/python_3.10")
  skip_if(python_path == "", "Python not found in system path")

  # 执行函数，捕获任何错误以进行测试
  result <- tryCatch({
    run_deconvolution_methods(scrna_path, spatial_path, 'celltype', python_path, output_path)
    TRUE
  }, error = function(e) {
    message("Error in run_deconvolution_methods: ", e$message)
    FALSE
  })

  # 验证函数是否成功运行
  expect_true(result, "run_deconvolution_methods completed without errors")

  # 检查输出文件是否创建
  expected_files <- c(
    file.path(output_path, "CARD_result.Rdata"),
    file.path(output_path, "RCTD_result.Rdata"),
    file.path(output_path, "SPOTlight_result.Rdata"),
    file.path(output_path, "SpatialDWLS_result.Rdata"),
    file.path(output_path, "AdRoit_result.Rdata"),
    file.path(output_path, "Redeconve_result.Rdata"),
    file.path(output_path, "Seurat_result.Rdata"),
    file.path(output_path, "SpatialDecon_result.Rdata")
  )

  for (file in expected_files) {
    expect_true(file.exists(file), paste("Output file not created:", file))
  }
})
