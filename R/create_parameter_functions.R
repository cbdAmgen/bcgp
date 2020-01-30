#' Create a list with random parameter values.
#'
#' \code{create_parameter_list} returns a list that contains randomly generated
#' parameter values.
#'
#' This creates a list that contains randomly generated parameter values for a
#' model corresponding to inputs \code{composite}, \code{stationary},
#' \code{noise}, and \code{d}. This function is meant to set up the parameters
#' before a call to \code{simulate_from_model()} or to \code{bcgpsims()}. The
#' user can modify the list of parameter values if desired. If a non-stationary
#' model is desired, the parameters that define the variance process,
#' \eqn{[\mu_V,\rho_V,\sigma^2_V]^\top}, will be returned and can then be
#' modified, but the variance process itself will not be returned and cannot be
#' specified.
#'
#' @param composite A logical, \code{TRUE} for a composite of a global process
#' and a local process, \code{FALSE} for non-composite. Defaults to \code{TRUE}.
#' @param stationary A logical, \code{FALSE} for a non-stationary process,
#' \code{TRUE} for a stationary process. If \code{FALSE}, the variance for the
#' process is \eqn{\sigma^2(x)}, and if \code{TRUE}, the variance is
#' \eqn{\sigma^2}. Defaults to \code{FALSE}.
#' @param noise If the data should be noise-free (such as from a deterministic
#' computer model), then \code{noise} should be \code{FALSE}. Otherwise, it
#' should be \code{TRUE}. Defaults to \code{FALSE}
#' @param d An integer giving the dimension of the data.
#' @return A list of randomly-generated parameter values corresponding to the
#' BCGP model described by \code{composite}, \code{stationary}, \code{noise},
#' and \code{d}.
#' @family preprocessing functions
#' @seealso \code{\link{simulate_from_model}} \code{\link{bcgpsims}}
#' @examples
#' create_parameter_list(composite = FALSE, stationary = TRUE, noise = FALSE,
#'                     d = 2)
#' create_parameter_list(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                     d = 1)
#' @export
create_parameter_list <- function(composite = TRUE, stationary = FALSE,
                                  noise = FALSE, d = 1L){

  d <- as.integer(d)

  if(isTRUE(composite)){
    if(isFALSE(stationary)){
      param_list <- create_param_comp_ns(d)
    }else{ # composite == TRUE, stationary == TRUE
      param_list <- create_param_comp_s(d)
    }
  }else{
    if(isFALSE(stationary)){
      param_list <- create_param_noncomp_ns(d)
    }else{ # composite == FALSE, stationary == TRUE
      param_list <- create_param_noncomp_s(d)
    }
  }
  param_list$sigma2eps <- ifelse(isTRUE(noise), rgamma(1, 2, scale = 0.01), 0)
  return(param_list)
}

create_param_comp_ns <- function(d){

  rhoG <- runif(d)
  rhoL <- runif(d, 0, rhoG)
  muV <- -0.1
  sigma2V <- 0.1
  rhoV <- runif(d)

  param_list <- list(beta0 = 0,
                     w = runif(1, 0.5, 1),
                     rhoG = rhoG,
                     rhoL = rhoL,
                     muV = muV,
                     sigma2V = sigma2V,
                     rhoV = rhoV)
  return(param_list)

}

create_param_comp_s <- function(d){

  rhoG <- runif(d)
  rhoL <- runif(d, 0, rhoG)

  param_list <- list(beta0 = 0,
                     w = runif(1, 0.5, 1),
                     rhoG = rhoG,
                     rhoL = rhoL,
                     sigma2 = 1)
  return(param_list)
}

create_param_noncomp_ns <- function(d){

  muV <- -0.1
  sigma2V <- 0.1
  rhoV <- runif(d)

  param_list <- list(beta0 = 0,
                     rho = runif(d),
                     muV = muV,
                     sigma2V = sigma2V,
                     rhoV = rhoV)
  return(param_list)
}

create_param_noncomp_s <- function(d){

  param_list <- list(beta0 = 0,
                     rho = runif(d),
                     sigma2 = 1)
  return(param_list)
}
