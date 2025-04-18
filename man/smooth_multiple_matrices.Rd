% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networksmooth.r
\name{smooth_multiple_matrices}
\alias{smooth_multiple_matrices}
\title{Smooth Multiple Matrices Using Autoencoder}
\usage{
smooth_multiple_matrices(
  proportion_matrices,
  pos,
  n_pos_features = 16,
  epochs = 4000,
  batch_size = 5,
  cropping = list(c(0, 9), c(9, 1))
)
}
\arguments{
\item{proportion_matrices}{A list of proportion matrices to be smoothed.}

\item{pos}{A data frame with x and y coordinates for each pixel.}

\item{n_pos_features}{Number of positional encoding features. Default is 16.}

\item{epochs}{Number of training epochs. Default is 4000.}

\item{batch_size}{Batch size for training. Default is 5.}

\item{cropping}{A list specifying the cropping parameters. Default is list(c(0, 9), c(9, 1)).}
}
\value{
A list containing the trained model, training history, and smoothed matrices.
}
\description{
This function smooths multiple proportion matrices using an autoencoder model.
}
