#' Create a correlation matrix.
#'
#' \code{get_cor_mat_R} returns a correlation matrix.
#'
#' This creates a correlation matrix, \code{R}, where \deqn{R_{ij} = \prod_{k =
#' 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}
#'
#' @param x An \code{n x d} matrix, where \code{d} is the dimension of the data.
#' @param rho A vector of length d. Each value in \code{rho} should be between
#'  0 and 1.
#' @return An \code{n x n} correlation matrix
#' @family correlation and covariance functions
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rho <- runif(d, 0, 1)
#' get_cor_mat_R(x, rho)
#' @export

get_cor_mat_R <- function(x, rho){

  # If I choose not to export this function, then I'll skip the error-checking
  # since the only time this function would be called is if everything is
  # correct.
  stopifnot(is.matrix(x), is.numeric(x), all(0 <= rho & rho <= 1),
            ncol(x) == length(rho))

  n <- nrow(x)
  d <- ncol(x)

  R <- matrix(0, nrow = n, ncol = n)

  for(i in 1:(n-1)){
    R[i, i] <- 1.0
    for(j in (i+1):n){
      R[i, j] <- 1.0
      dist <- x[i, ] - x[j, ]
      R[i, j] <- prod(rho^(16*dist^2))
      R[j, i] <- R[i, j]
    }

  }
  R[n, n] <- 1.0
  return(R)
}


#' Create a covariance matrix.
#'
#' \code{get_cov_mat_ns_R} returns a non-stationary covariance matrix (with
#' nugget).
#'
#' This creates a covariance matrix, \code{C}, where
#' \deqn{C =  V^{0.5}RV^{0.5} + \sigma^2_\epsilon I} where V is a matrix with
#' process variances (from vector V above) on the diagonal,
#' \code{R} is a correlation matrix, and \eqn{\sigma^2_\epsilon} is the variance
#'  of the noise (or nugget).
#'
#' @param V A positive vector of length \code{n}.
#' @param R An \code{n x n} correlation matrix.
#' @param sigma2eps A positive scalar representing the variance of the noise
#' (or nugget).
#'
#' @section Note:
#' Surprisingly, this method is substantially faster, even for relatively large
#' matrices (tested up to \code{5000 x 5000}), than doing sparse matrix
#' multiplication in the Matrix package. Sparse matrix multiplication was orders
#' of magnitude slower for small matrices than the method implemented above or
#' for \code{diag(V)^0.5 \%*\% R \%*\% diag(V)^0.5}, which was roughly the same
#' speed as above for small matrices, but much slower for large matrices.
#'
#' @return An \code{n x n} covariance matrix
#' @family correlation and covariance functions
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rho <- runif(d, 0, 1)
#' R <- get_cor_mat_R(x, rho)
#' sigma2eps <- 0.01
#' V <- rlnorm(n, -0.1, 0.1)
#' get_cov_mat_ns_R(V, R, sigma2eps)
#' @export
get_cov_mat_ns_R <- function(V, R, sigma2eps){

  # If I choose not to export this function, then I'll skip the error-checking
  # since the only time this function would be called is if everything is
  # correct.
  # stopifnot(all(V > 0),
  #           check_valid_cor_mat(R),
  #           (sig2eps >= 0))

  rootV <- sqrt(V)
  C <- t(rootV*R) * rootV + diag(sigma2eps, length(V))
  return(C)

  # NOTE: Surprisingly, this method is substantially faster, even for relatively
  # large matrices (tested up to 5000 x 5000), than doing sparse matrix
  # multiplication in the Matrix package.

  # Sparse matrix multiplication was orders of magnitude slower for small
  # matrices than the method implemented above or for
  # diag(V)^0.5 %*% R %*% diag(V)^0.5, which was roughly the same speed as above
  # for small matrices, but much slower for large matrices.

}

#' Create a covariance matrix.
#'
#' \code{get_cov_mat_s_R} returns a stationary covariance matrix (with nugget).
#' This creates a covariance matrix, C, where \deqn{C =  sigma^2 R + \sig2Eps I}
#' \code{R} is a correlation matrix, \code{sigma2} the process variance, and
#' \code{sigma2eps} is the variance of the noise (or nugget).
#' @param sigma2 A positive scalar representing the process variance
#' @param R An \emph{n x n} correlation matrix.
#' @param sigma2eps A positive scalar representing the variance of the noise
#' (or nugget).
#' @return An \code{n x n} covariance matrix
#' @family correlation and covariance functions
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rho <- runif(d, 0, 1)
#' R <- get_cor_mat_R(x, rho)
#' sigma2 <- 1.2
#' get_cov_mat_s_R(sigma2, R, sigma2eps)
#' @export
get_cov_mat_s_R <- function(sigma2, R, sigma2eps){
  sigma2 * R + diag(sigma2eps, nrow = nrow(R))
}

#' Combine correlation matrices.
#'
#' \code{combine_cor_mats_R} returns a correlation matrix that is a weighted sum
#'  of two other correlation matrices.
#'
#' This creates a correlation matrix, \emph{R}, where \deqn{R = wG + (1-w)L}. In
#' the \code{bcgp} setting, \emph{G} is the global correlation matrix, \emph{L}
#' is the local correlation matrix, and \emph{w} is the weight.
#'
#' @param w A scalar between 0 and 1. In the \code{bcgp} setting, the user will
#' have set a lower and upper bound (the default is to be between 0.5 and 1).
#' @param G An \emph{n x n} correlation matrix, often the result of
#' \code{\link{get_cor_mat_R}}.
#' @param L An \emph{n x n} correlation matrix, often the result of
#' \code{\link{get_cor_mat_R}}.
#' @return An \emph{n x n} correlation matrix
#' @family correlation and covariance functions
#' @examples
#' n <- 10
#' d <- 2
#' x <- matrix(runif(n * d), nrow = n, ncol = d)
#' rhoG <- runif(d, 0, 1)
#' rhoL <- runif(d, 0, rhoG)
#' G <- get_cor_mat_R(x, rhoG)
#' L <- get_cor_mat_R(x, rhoL)
#' w <- runif(1, 0.5, 1)
#' combine_cor_mats_R(w, G, L)
#' @export

combine_cor_mats_R <- function(w, G, L){

  # If I choose not to export this function, then I'll skip the error-checking
  # since the only time this function would be called is if everything is
  # correct.
  # stopifnot((0 <= w && w <= 1),
  #           check_valid_cor_mat(G),
  #           check_valid_cor_mat(L))

  R <- w*G + (1 - w)*L
  return(R)
}
