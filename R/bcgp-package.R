#' Bayesian Composite Gaussian Process Modeling
#'
#' @description The \pkg{bcgp} package enables Gaussian process (GP) regression
#' models to be estimated using Markov Chain Monte Carlo (MCMC).
#'
#' @details The \pkg{bcgp} package supports both stationary and non-stationary
#' models. For non-stationary models, the non-stationarity is introduced by
#' having a non-constant variance. This package also supports composite and
#' non-composite models, where \emph{composite} models are a sum of a global GP
#' that captures the overall trend and a local GP that makes small adjustments
#' around the trend.
#'
#' @docType package
#' @name bcgp-package
#' @aliases bcgp
#' @useDynLib bcgp, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Prediction Using a Bayesian Heteroscedastic Composite Gaussian Process (2020)
#'
NULL
