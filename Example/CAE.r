setwd("D:/Bioinfomation_Project/Create_R_Package")

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

py_config()

load("data/TestData/seqFISH+/CARD_result.Rdata")
CARD_matrix <- as.matrix(CARD_result)
load("data/TestData/seqFISH+/RCTD_result.Rdata")
RCTD_matrix <- as.matrix(RCTD_result)
load("data/TestData/seqFISH+/SPOTlight_result.Rdata")
SPOTlight_matrix <- as.matrix(SPOTlight_result)
load("data/TestData/seqFISH+/SpatialDWLS_result.Rdata")
SpatialDWLS_matrix <- as.matrix(SpatialDWLS_result)
load(file = 'data/TestData/seqFISH+/pos.Rdata')
pos2 <- pos

# Make the column names of the matrix names the same
colnames(SPOTlight_matrix) = gsub("\\.", "_", colnames(SPOTlight_matrix))

CARD_matrix <- CARD_matrix[,colnames(RCTD_matrix)]
SPOTlight_matrix <- SPOTlight_matrix[,colnames(RCTD_matrix)]
SpatialDWLS_matrix <- SpatialDWLS_matrix[,colnames(RCTD_matrix)]


row_names <- rownames(SpatialDWLS_matrix)
col_names <- colnames(SpatialDWLS_matrix)
SpatialDWLS_matrix <- matrix(as.numeric(SpatialDWLS_matrix), nrow = nrow(SpatialDWLS_matrix), ncol = ncol(SpatialDWLS_matrix))
rownames(SpatialDWLS_matrix) <- row_names
colnames(SpatialDWLS_matrix) <- col_names
CARD_matrix
RCTD_matrix
SPOTlight_matrix
SpatialDWLS_matrix

# Prepare data
proportion_matrices <- list(CARD_matrix, RCTD_matrix, SPOTlight_matrix, SpatialDWLS_matrix)

result <- smooth_multiple_matrices(
  proportion_matrices, 
  pos2,
  epochs = 4000,
  batch_size = 5,
  cropping = list(c(0, 9), c(9, 1))  # Custom Cropping Parameters
)

# Get trained models and training history
trained_model <- result$model
training_history <- result$history
# smoothed_matrices <- result$smoothed_matrices
smoothed_matrix <- result$smoothed_matrices
plot(training_history)


smoothed_matrix
save(smoothed_matrix, file = "data/TestData/seqFISH+/smoothed_matrix.RData")
