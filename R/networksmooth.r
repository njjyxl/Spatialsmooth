# Sys.setenv(RETICULATE_PYTHON = "D:/Application/anaconda3/envs/python_3.10")

# options(warn = -1)
# suppressMessages(library(spacexr))
# suppressMessages(library(Matrix))
# suppressMessages(library(data.table))
# suppressMessages(library(Seurat))
# suppressMessages(library(SeuratDisk))
# suppressMessages(library(gstat))
# suppressMessages(library(sp))
# suppressMessages(library(igraph))
# suppressMessages(library(fields))
# suppressMessages(library(keras))
# suppressMessages(library(tensorflow))
# suppressMessages(library(reticulate))

# py_config()

#' Smooth Multiple Matrices Using Autoencoder
#'
#' This function smooths multiple proportion matrices using an autoencoder model.
#'
#' @param proportion_matrices A list of proportion matrices to be smoothed.
#' @param pos A data frame with x and y coordinates for each pixel.
#' @param n_pos_features Number of positional encoding features. Default is 16.
#' @param epochs Number of training epochs. Default is 4000.
#' @param batch_size Batch size for training. Default is 5.
#' @param cropping A list specifying the cropping parameters. Default is list(c(0, 9), c(9, 1)).
#'
#' @return A list containing the trained model, training history, and smoothed matrices.
#' @export
#'
#' @importFrom keras layer_input layer_conv_2d layer_max_pooling_2d layer_dropout layer_upsampling_2d layer_cropping_2d keras_model
#' @importFrom tensorflow k_abs k_mean k_square
smooth_multiple_matrices <- function(proportion_matrices, pos, n_pos_features = 16, epochs = 4000, batch_size = 5, cropping = list(c(0, 9), c(9, 1))) {
  # 确保所有矩阵有相同的维度
  n_samples <- length(proportion_matrices)
  n_pixels <- nrow(proportion_matrices[[1]])
  n_cell_types <- ncol(proportion_matrices[[1]])
  
  # 位置编码函数
  pos_encoding <- function(pos, freq) {
    angle_rads_x <- matrix(pos$x, ncol = length(freq), nrow = nrow(pos)) * freq
    angle_rads_y <- matrix(pos$y, ncol = length(freq), nrow = nrow(pos)) * freq
    
    angle_rads_x[, seq(2, length(freq), 2)] <- cos(angle_rads_x[, seq(2, length(freq), 2)])
    angle_rads_x[, seq(1, length(freq), 2)] <- sin(angle_rads_x[, seq(1, length(freq), 2)])
    
    angle_rads_y[, seq(2, length(freq), 2)] <- cos(angle_rads_y[, seq(2, length(freq), 2)])
    angle_rads_y[, seq(1, length(freq), 2)] <- sin(angle_rads_y[, seq(1, length(freq), 2)])
    
    pos_encoding <- cbind(angle_rads_x, angle_rads_y)
    
    return(pos_encoding)
  }
  
  # 执行位置编码
  freq <- 2^seq(0, (n_pos_features/2 - 1))
  pos_encoded <- pos_encoding(pos, freq)
  
  n_features <- ncol(pos_encoded) + 1
  x_train <- array(0, dim = c(n_samples, n_pixels, n_cell_types, n_features))
  
  for (i in 1:n_samples) {
    x_train[i, , , 1] <- proportion_matrices[[i]]
    for (j in 1:n_pixels) {
      x_train[i, j, , 2:n_features] <- matrix(rep(pos_encoded[j, ], n_cell_types), ncol = ncol(pos_encoded), byrow = TRUE)
    }
  }
  
  # 定义模型架构
  input_img <- layer_input(shape = dim(x_train)[-1])
  
  # Encoder
  x <- input_img %>% 
    layer_conv_2d(filters = 32, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_max_pooling_2d(pool_size = 2, padding = "same") %>%
    layer_dropout(rate = 0.25) %>%
    layer_conv_2d(filters = 64, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_max_pooling_2d(pool_size = 2, padding = "same") %>%
    layer_dropout(rate = 0.25) %>%
    layer_conv_2d(filters = 128, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_max_pooling_2d(pool_size = 2, padding = "same") %>%
    layer_dropout(rate = 0.25) %>%
    layer_conv_2d(filters = 256, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_max_pooling_2d(pool_size = 2, padding = "same")
  
  # Bottleneck
  x <- x %>%
    layer_conv_2d(filters = 512, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_dropout(rate = 0.5)
  
  # Decoder
  x <- x %>%
    layer_conv_2d(filters = 256, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_upsampling_2d(size = 2) %>%
    layer_conv_2d(filters = 128, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_upsampling_2d(size = 2) %>%
    layer_dropout(rate = 0.25) %>%
    layer_conv_2d(filters = 64, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_upsampling_2d(size = 2) %>%
    layer_dropout(rate = 0.25) %>%
    layer_conv_2d(filters = 32, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_upsampling_2d(size = 2) %>%
    layer_dropout(rate = 0.25) %>%
    layer_conv_2d(filters = 16, kernel_size = 3, padding = "same", activation = "relu") %>%
    layer_conv_2d(filters = n_features, kernel_size = 3, padding = "same", activation = "softmax") %>%
    layer_cropping_2d(cropping = cropping)  # 使用传入的cropping参数
  
  # 定义自动编码器
  autoencoder <- keras_model(input_img, x)
  
  # 定义损失函数
  total_variation_loss <- function(x) {
    h_diff <- k_abs(x[, 1:(n_pixels - 1), , ] - x[, 2:n_pixels, , ])
    v_diff <- k_abs(x[, , 1:(n_cell_types - 1), ] - x[, , 2:n_cell_types, ])
    tv_loss <- k_mean(h_diff) + k_mean(v_diff)
    return(tv_loss)
  }
  
#   multi_task_loss <- function(y_true, y_pred) {
#     mse_loss <- k_mean(k_square(y_true[,,, 1] - y_pred[,,, 1]))
#     tv_loss <- total_variation_loss(y_pred)
#     lambda_mse <- 1.0
#     lambda <- 0.001  
#     pos_loss <- k_mean(k_square(y_true[,,, 2:n_features] - y_pred[,,, 2:n_features]))
#     total_loss <- lambda_mse * mse_loss + lambda * tv_loss + 0.001 * pos_loss
#     return(total_loss)
#   }
  weighted_mse_with_tv <- function(y_true, y_pred) {
    # Assume the weight matrix is `weights`, with the same spatial dimensions as the input data
    weights <- k_ones_like(y_true)
#   weights <- create_weight_tensor(y_true)
    # Calculate weighted mean squared error
    squared_difference <- k_square(y_true - y_pred)
    weighted_squared_difference <- k_sum(squared_difference * weights)
    wmse <- weighted_squared_difference / k_sum(weights)
  
    # Calculate total variation regularization term
    tv_loss <- total_variation_loss(y_pred)
    # 动态权重
#    lambda <- k_exp(-k_get_value(autoencoder$optimizer$iterations) / 1000) * 0.1
    
#   pos_loss <- k_mean(k_square(y_true[,,, 2:n_features] - y_pred[,,, 2:n_features]))
    lambda <- 0.01
    # Add weighted mean squared error and total variation regularization term, and convert to scalar value
#   total_loss <- k_mean(wmse + lambda * tv_loss + 0.001 * pos_loss)
    total_loss <- k_mean(wmse + lambda * tv_loss)
    return(total_loss)
  }
  
  # 编译模型
#   autoencoder %>% compile(optimizer = 'adam', loss = multi_task_loss)
  autoencoder %>% compile(optimizer = 'adam', loss = weighted_mse_with_tv)
  
  # 定义回调
  early_stopping <- callback_early_stopping(
    monitor = "val_loss", 
    patience = 100, 
    restore_best_weights = TRUE
  )
  
  checkpoint <- callback_model_checkpoint(
    filepath = "best_model.h5",
    monitor = "val_loss",
    save_best_only = TRUE
  )
  
  # 训练模型
  history <- autoencoder %>% fit(
    x_train, x_train,
    epochs = epochs,
    batch_size = batch_size,
    shuffle = TRUE,
    validation_split = 0.2,
    callbacks = list(early_stopping, checkpoint)
  )
  # Predict the smoothed proportions
  smoothed_prop <- autoencoder %>% predict(x_train)
  
  # Reshape the smoothed proportions back to the original shape
  smoothed_prop <- array_reshape(smoothed_prop, dim = c(n_samples, n_pixels, n_cell_types, n_features))

  # Extract the smoothed proportion matrices
  smoothed_matrices <- lapply(1:n_samples, function(i) smoothed_prop[i, , , 1])

  # Normalize the smoothed proportions to sum up to 1 for each pixel
  smoothed_matrices <- lapply(smoothed_matrices, function(mat) apply(mat, 1, function(x) x / sum(x)))

  matrices_transposed <- lapply(smoothed_matrices, t)

  # Assign row and column names to the smoothed matrices
  for (i in 1:n_samples) {
    rownames(matrices_transposed[[i]]) <- rownames(proportion_matrices[[i]])
    colnames(matrices_transposed[[i]]) <- colnames(proportion_matrices[[i]])
  }
  ssim <- function(x, y) {
  mu_x <- mean(x)
  mu_y <- mean(y)
  var_x <- var(x)
  var_y <- var(y)
  covar_xy <- cov(c(x), c(y))
  
  c1 <- (0.01 * 255)^2
  c2 <- (0.03 * 255)^2
  
  num <- (2 * mu_x * mu_y + c1) * (2 * covar_xy + c2)
  den <- (mu_x^2 + mu_y^2 + c1) * (var_x + var_y + c2)
  
  return(num / den)
  }

  calculate_tv <- function(matrix) {
    h_diff <- abs(matrix[, 1:(ncol(matrix)-1)] - matrix[, 2:ncol(matrix)])
    v_diff <- abs(matrix[1:(nrow(matrix)-1), ] - matrix[2:nrow(matrix), ])
    sum(h_diff) + sum(v_diff)
  }

  scores <- sapply(1:length(proportion_matrices), function(i) {
    smoothed <- matrices_transposed[[i]]
    original <- proportion_matrices[[i]]
    cor_value <- cor(as.vector(smoothed), as.vector(original))
    ssim_value <- ssim(as.vector(smoothed), as.vector(original))
    tv_value <- calculate_tv(original) - calculate_tv(smoothed)
    c(cor = cor_value, ssim = ssim_value, TV = tv_value)
  })

  cor_values <- scores["cor", ]
  ssim_values <- scores["ssim", ]
  tv_values <- scores["TV", ]

  # Normalize the scores
  normalized_cor <- (cor_values - min(cor_values)) / (max(cor_values) - min(cor_values))
  normalized_ssim <- (ssim_values - min(ssim_values)) / (max(ssim_values) - min(ssim_values))

  # Calculate combined score (equal weight for both metrics)
  combined_scores <- (normalized_cor + normalized_ssim + 0.05 * tv_values) / 3

  # Find the index of the matrix with the highest combined score
  best_index <- which.max(combined_scores)

  # Select the best matrix
  best_matrix <- matrices_transposed[[best_index]]

  # Assign row and column names to the best matrix
  rownames(best_matrix) <- rownames(proportion_matrices[[best_index]])
  colnames(best_matrix) <- colnames(proportion_matrices[[best_index]])

  # Print the best matrix
  cat("Best Smoothed Matrix (from method", best_index, "):\n")
  print(best_matrix)
  smoothed_matrix <- best_matrix

  # Print which method produced the best result
  best_method <- c("CARD", "RCTD", "SPOTlight", "SpatialDWLS")[best_index]
  # best_method <- c("CARD", "RCTD")[best_index]
  cat("Best method:", best_method, "\n")

#   # Print the scores
#   cat("Correlation coefficient:", cor_values[best_index], "\n")
#   cat("SSIM:", ssim_values[best_index], "\n")
#   cat("TV:", tv_values[best_index], "\n")
  cat("Combined score:", combined_scores[best_index], "\n")

  # Return the results
  return(list(
    model = autoencoder,
    history = history,
    smoothed_matrices = smoothed_matrix
  ))
  # 返回训练好的模型和训练历史
#   return(list(model = autoencoder, history = history))
}