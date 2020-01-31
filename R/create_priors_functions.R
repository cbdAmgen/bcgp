#' Create a list of priors for BCGP
#'
#' \code{create_priors} returns a list of default priors.
#'
#' This creates a list that contains default hyperparameter values for a
#' model corresponding to inputs \code{composite}, \code{stationary},
#' \code{noise}, and \code{d}. This function is meant to set up the priors in a
#' \code{bcgpmodel} object before a call to \code{bcgp_sampling()}. The
#' user will be able to modify the hyperparameter values in the \code{bcgpmodel}
#' object before the call to \code{bcgp_sampling()} if desired, but the user
#' has no choice of prior distributions (only the parameters that define the
#' given distribution).
#'
#' @param composite A logical, \code{TRUE} for a composite of a global process,
#' a local process, and an error process, \code{FALSE} for non-composite.
#' Defaults to \code{TRUE}.
#' @param stationary A logical, \code{FALSE} for a non-stationary process,
#' \code{TRUE} for a stationary process. If \code{FALSE}, the variance for the
#' process is \eqn{\sigma^2(x)}, and if \code{TRUE}, the variance is
#' \eqn{\sigma^2}. Defaults to \code{FALSE}.
#' @param noise If the data should be noise-free (such as from a deterministic
#' computer model), then \code{noise} should be \code{FALSE}. Otherwise, it
#' should be \code{TRUE}. Defaults to \code{FALSE}
#' @param d An integer giving the dimension of the data.
#' @return A list containing the default values for all the prior parameters,
#' information about the process, and information about the prior distributions.
#' @seealso \linkS4class{bcgpmodel} \code{\link{bcgpfit}}
#' \code{\link{bcgp}}
#' @examples
#' create_priors(composite = TRUE, stationary = FALSE, noise = FALSE, d = 1)
#' create_priors(composite = FALSE, stationary = TRUE, noise = TRUE, d = 3)

create_priors <- function(composite = TRUE, stationary = FALSE,
                          noise = FALSE, d = 1L){

  d <- as.integer(d)

  if(isTRUE(composite)){
    if(isFALSE(stationary)){
      ## composite, non- stationary
      prior_info <- create_priors_comp_ns(d)
    }else{
      ## composite, stationary
      prior_info <- create_priors_comp_s(d)
    }
  }else{
    if(isFALSE(stationary)){
      ## non-composite, non- stationary
      prior_info <- create_priors_noncomp_ns(d)
    }else{
      ## non-composite, stationary
      prior_info <- create_priors_noncomp_s(d)
    }
  }

  if(isFALSE(noise)){
    prior_info$priors$sigma2eps <- list(alpha = 1.1,
                                        beta = 1e-6)
  }else{
    prior_info$priors$sigma2eps <- list(alpha = 1e-1,
                                        beta = 1e-1)
  }

  return(prior_info)
}

create_priors_comp_ns <- function(d){
  prior_list <- list(w = list(lower = 0.5,
                              upper = 1.0,
                              alpha = 1,
                              beta = 1),
                     rhoG = list(alpha = rep(1, d),
                                 beta = rep(1, d)),
                     rhoL = list(alpha = rep(1, d),
                                 beta = rep(1, d)),
                     muV = list(betaV = -0.1,
                                sigma2 = 0.1),
                     rhoV = list(alpha = rep(1, d),
                                 beta = rep(1, d)),
                     sigma2V = list(alpha = 2 + sqrt(0.1),
                                    beta = 100/(1+sqrt(1/10))))
  prior_dists <- list(beta0 = "noninformative",
                      w = "TrBeta(lower, upper, alpha, beta)",
                      rhoG = "Beta(alpha, beta)",
                      rhoL = "TrBeta(0, rhoG, alpha, beta)",
                      sigma2eps = "Gamma(alpha, scale = beta)",
                      muV = "Lognormal(betaV, variance = sigma2)",
                      rhoV = "Beta(alpha, beta)",
                      sigma2V = "Inverse Gamma(alpha, beta)")

  return(list(priors = prior_list, distributions = prior_dists))
}

create_priors_comp_s <- function(d){
  prior_list <- list(w = list(lower = 0.5,
                              upper = 1.0,
                              alpha = 1,
                              beta = 1),
                     rhoG = list(alpha = rep(1, d),
                                 beta = rep(1, d)),
                     rhoL = list(alpha = rep(1, d),
                                 beta = rep(1, d)),
                     sigma2 = list(alpha = 1,
                                   beta = 1))
  prior_dists <- list(beta0 = "noninformative",
                      w = "TrBeta(lower, upper, alpha, beta)",
                      rhoG = "Beta(alpha, beta)",
                      rhoL = "TrBeta(0, rhoG, alpha, beta)",
                      sigma2eps = "Gamma(alpha, scale = beta)",
                      sigma2 = "Gamma(alpha, scale = beta)")

  return(list(priors = prior_list, distributions = prior_dists))
}

create_priors_noncomp_ns <- function(d){
  prior_list <- list(rho = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                     muV = list(betaV = -0.1,
                                sigma2 = 0.1),
                     rhoV = list(alpha = rep(1, d),
                                 beta = rep(1, d)),
                     sigma2V = list(alpha = 2 + sqrt(0.1),
                                    beta = 100/(1+sqrt(1/10))))
  prior_dists <- list(beta0 = "noninformative",
                      rho = "Beta(alpha, beta)",
                      sigma2eps = "Gamma(alpha, beta)",
                      muV = "Lognormal(betaV, variance = sigma2)",
                      rhoV = "Beta(alpha, beta)",
                      sigma2V = "Inverse Gamma(alpha, beta)")

  return(list(priors = prior_list, distributions = prior_dists))
}

create_priors_noncomp_s <- function(d){
  prior_list <- list(rho = list(alpha = rep(1, d),
                                beta = rep(1, d)),
                     sigma2 = list(alpha = 1,
                                   beta = 1))
  prior_dists <- list(beta0 = "noninformative",
                      rho = "Beta(alpha, beta)",
                      sigma2eps = "Gamma(alpha, beta)",
                      sigma2 = "Gamma(alpha, beta)")

  return(list(priors = prior_list, distributions = prior_dists))
}
