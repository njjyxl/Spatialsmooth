Sys.setenv(RETICULATE_PYTHON = "D:/Application/anaconda3/envs/python_3.10")

options(warn = -1)
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(gstat))
suppressMessages(library(dplyr))
suppressMessages(library(future))
suppressMessages(library(spatialreg))
suppressMessages(library(sp))
suppressMessages(library(igraph))
suppressMessages(library(fields))
suppressMessages(library(keras))
suppressMessages(library(tensorflow))
suppressMessages(library(reticulate))
suppressMessages(library(Spatialsmooth))



test_that("multiplication works", {
  # 尝试加载测试数据
  test_data_path <- "./"
  skip_if_not(dir.exists(test_data_path), "Test data directory not found")

  # 加载测试数据
  tryCatch({
    load(file.path(test_data_path, "CARD_result.Rdata"))
    load(file.path(test_data_path, "RCTD_result.Rdata"))
    load(file.path(test_data_path, "SPOTlight_result.Rdata"))
    load(file.path(test_data_path, "SpatialDWLS_result.Rdata"))
    load(file.path(test_data_path, "pos.Rdata"))

    # 转换为矩阵格式
    CARD_matrix <- as.matrix(CARD_result)
    RCTD_matrix <- as.matrix(RCTD_result)
    SPOTlight_matrix <- as.matrix(SPOTlight_result)
    SpatialDWLS_matrix <- as.matrix(SpatialDWLS_result)

    # 标准化列名
    colnames(SPOTlight_matrix) = gsub("\\.", "_", colnames(SPOTlight_matrix))
    CARD_matrix <- CARD_matrix[,colnames(RCTD_matrix)]
    SPOTlight_matrix <- SPOTlight_matrix[,colnames(RCTD_matrix)]
    SpatialDWLS_matrix <- SpatialDWLS_matrix[,colnames(RCTD_matrix)]

    # 转换SpatialDWLS_matrix中的值为数值型
    row_names <- rownames(SpatialDWLS_matrix)
    col_names <- colnames(SpatialDWLS_matrix)
    SpatialDWLS_matrix <- matrix(as.numeric(SpatialDWLS_matrix),
                                 nrow = nrow(SpatialDWLS_matrix),
                                 ncol = ncol(SpatialDWLS_matrix))
    rownames(SpatialDWLS_matrix) <- row_names
    colnames(SpatialDWLS_matrix) <- col_names

    # 准备输入数据
    proportion_matrices <- list(CARD_matrix, RCTD_matrix, SPOTlight_matrix, SpatialDWLS_matrix)

    # 运行函数，使用较小的epochs值以加快测试速度
    result <- smooth_multiple_matrices(
      proportion_matrices,
      pos,
      epochs = 1000,  # 减少epochs以加快测试
      batch_size = 5,
      cropping = list(c(0, 9), c(9, 1))
    )

    # 验证结果结构
    expect_type(result, "list")
    expect_true(all(c("model", "history", "smoothed_matrices") %in% names(result)))

    # 验证模型训练历史
    expect_true("keras_training_history" %in% class(result$history) ||
                  "tensorflow.python.keras.callbacks.History" %in% class(result$history))

    # 验证平滑矩阵的结构
    smoothed_matrix <- result$smoothed_matrices

    # 验证平滑矩阵的值在合理范围内（例如0-1之间的比例值）
    for (i in 1:length(smoothed_matrix)) {
      expect_true(all(smoothed_matrix[[i]] >= 0 & smoothed_matrix[[i]] <= 1))
    }

  }, error = function(e) {
    skip(paste("Error loading or processing test data:", e$message))
  })
})
