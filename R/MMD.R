#' Internal computation of biased Gaussian MMD
#'
#' Computes the biased squared Maximum Mean Discrepancy (MMD)
#' between two samples using a Gaussian kernel with a diagonal
#' covariance matrix.
#'
#' @param y Numeric matrix. Rows represent observations.
#' @param x Numeric matrix. Rows represent observations.
#' @param var Numeric matrix representing the kernel variance-covariance matrix.
#'
#' @return A numeric value containing the biased squared MMD.
#'
#' @keywords internal
#' 
#' @details
#' This implementation is based on the Gaussian-kernel Maximum Mean
#' Discrepancy (MMD) formulation described by Gretton et al. (2007)
#' and follows the corresponding calculation implemented in the
#' EasyMMD R package by Anthony Ebert.
#' The implementation was rewritten for internal use in
#' robustICP and does not reproduce the original EasyMMD
#' source code. The implementation was validated against the original
#' EasyMMD implementation for equivalent inputs.
#'
#' @references
#' Gretton, A., Borgwardt, K. M., Rasch, M., Schölkopf, B., & Smola,
#' A. J. (2007). A kernel method for the two-sample-problem.
#' In Advances in Neural Information Processing Systems, 19, 513--520.
#'
#' Ebert, A. EasyMMD: Allows users to compute the Maximum Mean
#' Discrepancy (MMD) between two samples. R package, version 0.5.4.
#' \url{https://github.com/AnthonyEbert/EasyMMD}
.mmd_gaussian <- function(y, x, var) {
  
  stopifnot(is.matrix(y))
  stopifnot(is.matrix(x))
  stopifnot(ncol(y) == ncol(x))
  stopifnot(is.matrix(var))
  stopifnot(nrow(var) == ncol(x))
  stopifnot(ncol(var) == ncol(x))
  
  var_diag <- diag(var)
  
  if (any(!is.finite(var_diag)) || any(var_diag <= 0)) {
    stop(
      "The diagonal elements of 'var' must be positive and finite."
    )
  }
  
  # Sum of Gaussian kernel values
  #
  # This follows the same pairwise Mahalanobis-distance
  # calculation used by EasyMMD.
  kernel_sum <- function(a, b) {
    
    total <- 0
    
    for (i in seq_len(nrow(a))) {
      
      # Differences between observation i of a
      # and all observations of b
      differences <- sweep(
        b,
        2,
        a[i, ],
        FUN = "-"
      )
      
      # Squared Mahalanobis distances using diagonal var
      mahal_sq <- rowSums(
        sweep(
          differences^2,
          2,
          var_diag,
          FUN = "/"
        )
      )
      
      # Gaussian kernel
      total <- total +
        sum(exp(-0.5 * mahal_sq))
    }
    
    total
  }
  
  n_x <- nrow(x)
  n_y <- nrow(y)
  
  term_xx <- kernel_sum(x, x) / n_x^2
  term_yy <- kernel_sum(y, y) / n_y^2
  term_xy <- kernel_sum(x, y) / (n_x * n_y)
  
  term_xx + term_yy - 2 * term_xy
}