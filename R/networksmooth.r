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
  # Ensure all matrix have the same dimensions
  n_samples <- length(proportion_matrices)
  n_pixels <- nrow(proportion_matrices[[1]])
  n_cell_types <- ncol(proportion_matrices[[1]])
  
  # posing encoding
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
    layer_cropping_2d(cropping = cropping)  # Use the incoming cropping parameter
  
  
  autoencoder <- keras_model(input_img, x)
  
  # loss function
  total_variation_loss <- function(x) {
    h_diff <- k_abs(x[, 1:(n_pixels - 1), , ] - x[, 2:n_pixels, , ])
    v_diff <- k_abs(x[, , 1:(n_cell_types - 1), ] - x[, , 2:n_cell_types, ])
    tv_loss <- k_mean(h_diff) + k_mean(v_diff)
    return(tv_loss)
  }
  structural_similarity_loss <- function(y_true, y_pred) {
  mu_x <- k_mean(y_true, axis = c(2, 3), keepdims = TRUE)
  mu_y <- k_mean(y_pred, axis = c(2, 3), keepdims = TRUE)
  
  sigma_x <- k_mean(k_square(y_true - mu_x), axis = c(2, 3), keepdims = TRUE)
  sigma_y <- k_mean(k_square(y_pred - mu_y), axis = c(2, 3), keepdims = TRUE)
  sigma_xy <- k_mean((y_true - mu_x) * (y_pred - mu_y), axis = c(2, 3), keepdims = TRUE)
  
  c1 <- 0.01^2
  c2 <- 0.03^2
  
  ssim <- (2 * mu_x * mu_y + c1) * (2 * sigma_xy + c2) / 
          ((k_square(mu_x) + k_square(mu_y) + c1) * (sigma_x + sigma_y + c2))
  
  return(1 - k_mean(ssim))
  }
  
  multi_task_loss <- function(y_true, y_pred) {
    ss_loss <- structural_similarity_loss(y_true, y_pred)
    tv_loss <- total_variation_loss(y_pred)
    lambda <- 0.01
    total_loss <- k_mean( 0.5 * ss_loss + lambda * tv_loss)
    return(total_loss)
  }
  
  
  autoencoder %>% compile(optimizer = 'adam', loss = multi_task_loss)
  
  early_stopping <- callback_early_stopping(
    monitor = "val_loss", 
    patience = 300, 
    restore_best_weights = TRUE
  )
  
  checkpoint <- callback_model_checkpoint(
    filepath = "best_model.h5",
    monitor = "val_loss",
    save_best_only = TRUE
  )
  
  # training model
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
  smoothed_prop <- array_reshape(smoothed_prop, dim = c(n_samples, n_pixels, n_cell_types, n_features))
  smoothed_matrices <- lapply(1:n_samples, function(i) smoothed_prop[i, , , 1])
  smoothed_matrices <- lapply(smoothed_matrices, function(mat) apply(mat, 1, function(x) x / sum(x)))
  matrices_transposed <- lapply(smoothed_matrices, t)

  
  for (i in 1:n_samples) {
    rownames(matrices_transposed[[i]]) <- rownames(proportion_matrices[[i]])
    colnames(matrices_transposed[[i]]) <- colnames(proportion_matrices[[i]])
  }
  ssim <- function(im1, im2, M = 1) {
   
   im1 <- im1 / max(im1)
   im2 <- im2 / max(im2)
  
   mu1 <- mean(im1)
   mu2 <- mean(im2)
  
   sigma1 <- sqrt(mean((im1 - mu1) ^ 2))
   sigma2 <- sqrt(mean((im2 - mu2) ^ 2))
  
   sigma12 <- mean((im1 - mu1) * (im2 - mu2))
  
   k1 <- 0.01
   k2 <- 0.03
   L <- M
   C1 <- (k1 * L) ^ 2
   C2 <- (k2 * L) ^ 2
   C3 <- C2 / 2
  
   l12 <- (2 * mu1 * mu2 + C1) / (mu1 ^ 2 + mu2 ^ 2 + C1)
   c12 <- (2 * sigma1 * sigma2 + C2) / (sigma1 ^ 2 + sigma2 ^ 2 + C2)
   s12 <- (sigma12 + C3) / (sigma1 * sigma2 + C3)
  
   ssim_value <- l12 * c12 * s12
   return(ssim_value)
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

  normalized_cor <- (cor_values - min(cor_values)) / (max(cor_values) - min(cor_values))
  normalized_ssim <- (ssim_values - min(ssim_values)) / (max(ssim_values) - min(ssim_values))

  combined_scores <- (1.5 * normalized_cor + normalized_ssim + 0.003 * tv_values) / 3
  best_index <- which.max(combined_scores)

  best_matrix <- matrices_transposed[[best_index]]
  rownames(best_matrix) <- rownames(proportion_matrices[[best_index]])
  colnames(best_matrix) <- colnames(proportion_matrices[[best_index]])

  cat("Best Smoothed Matrix (from method", best_index, "):\n")
  print(best_matrix)
  smoothed_matrix <- best_matrix

  best_method <- c("CARD", "RCTD", "SPOTlight", "SpatialDWLS")[best_index]
  
  cat("Best method:", best_method, "\n")
                                                                     
  cat("Correlation coefficient:", cor_values[best_index], "\n")
  cat("SSIM:", ssim_values[best_index], "\n")
  cat("TV:", tv_values[best_index], "\n")
  cat("Combined score:", combined_scores[best_index], "\n")

  # Return the results
  return(list(
    model = autoencoder,
    history = history,
    smoothed_matrices = smoothed_matrix
  ))
}