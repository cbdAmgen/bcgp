#' Simulate from the model
#'
#' This function creates an object of class bcgpsims, simulating data from the
#' specified (by stationary, composite, noise) Gaussian process model.
#'
#' \code{simulate_from_model} returns an instance of S4 class \code{bcgpsims}.
#' This object can then be plotted to get an idea of what draws from these
#' models look like, or the data in the object can be fit.
#' \code{bcgpsims()} can also be called to create a \code{bcgpsims}
#' object that contains simulated data
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
#' @param n An integer giving the desired number of training data locations.
#' @param n_test An integer giving the desired number of test data locations.
#' @param parameters A list containing desired parameter values. If missing,
#' then parameter values will be drawn at random. A call to
#' \code{create_parameter_list()} will assist in the correct creation of this
#' list.
#' @param seed A numeric value indicating the seed for random number generation.
#' \code{as.integer} will be applied to the value before setting the seed for
#' the random number generator. The default is generated from 1 to the maximum
#' integer supported by R on the machine.
#' @param ... optional parameters
#' \describe{
#'   \item{\code{randomX1D}}{A logical indicating whether the training data
#'   should be generated in a sequence, \code{seq(0, 1, length.out = n)}, or
#'   randomly generated from [0, 1]. Defaults to \code{FALSE} (sequence). Only
#'   useful for 1-D data.}
#'   \item{\code{grid_test}}{A logical indicating whether the test data should
#'   be generated on a grid. Defaults to \code{FALSE}. Only useful for
#'   \eqn{d \ge 2}.}
#'   \item{\code{grid_test_size}}{An integer indicating the number of points per
#'   dimension for the test grid. Only useful for \eqn{d \ge 2}.}
#'   \item{}{Be aware that a grid in high dimensions quickly gets very large.
#'   For example, for \code{d} = 3 and \code{grid_test_size} = 10, the number of
#'   points in this grid is \eqn{10^3 = 1000}. Therefore, in higher dimensions (
#'   \eqn{d > 4}), the simulation will default to \code{grid_test = FALSE} if
#'   \code{grid_test_size} is left unspecified, and \code{n_test} test locations
#'   will be randomly selected on \eqn{[0, 1]^d}.}
#' }
#' @return An instance of S4 class \code{bcgpsims}
#' @seealso \linkS4class{bcgpsims} \code{\link{bcgpsims}}
#' \code{\link{create_parameter_list}}
#' @examples
#' simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE)
#'
#' simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                     d = 2, n = 20, n_test = 80)
#'
#' params <- create_parameter_list()
#' params$w <- 0.99
#' simulate_from_model(parameters = params, randomX1D = TRUE)
#' @export
simulate_from_model <- function(composite = TRUE, stationary = FALSE,
                                noise = FALSE, d = 1L, n = 15*d, n_test = 100*d,
                                parameters = create_parameter_list(composite,
                                                                   stationary,
                                                                   noise, d),
                                seed = sample.int(.Machine$integer.max, 1),
                                ...){

  extra_args <- sims_unpack_dots(d, ...)

  seed <- check_seed(seed)
  set.seed(seed)

  x_matrices <- create_x_and_x_test(d, n, n_test,
                                    grid_test = extra_args$grid_test,
                                    randomX1D = extra_args$randomX1D,
                                    grid_test_size = extra_args$grid_test_size)
  x <- x_matrices$x
  x_test <- x_matrices$x_test
  rm(x_matrices)

  validate_parameter_list(parameters = parameters,
                          composite = composite,
                          stationary = stationary,
                          d = d)

  set.seed(seed)
  data <- simulateY(x = x, x_test = x_test, parameters = parameters,
                    stationary = stationary, composite = composite,
                    seed = seed)

  if(isFALSE(composite)){
    toReturn <- new("bcgpsims",
                    training = list(x = x, y = data$y),
                    test = list(x = x_test, y = data$y_test,
                                grid = extra_args$grid_test),
                    parameters = data$parameters,
                    stationary = stationary,
                    composite = composite,
                    seed = seed)
  }else{
    toReturn <- new("bcgpsims",
                    training = list(x = x, y = data$y, yG = data$yG,
                                    yL = data$yL, yE = data$yE),
                    test = list(x = x_test, y = data$y_test,
                                yG = data$yGTest, yL = data$yLTest,
                                yE = data$yETest, grid = extra_args$grid_test),
                    parameters = data$parameters,
                    stationary = stationary,
                    composite = composite,
                    seed = seed)
  }

  return(toReturn)
}

simulateY <- function(x, x_test, parameters,
                      stationary, composite,
                      decomposition, seed = seed){

  set.seed(seed)

  if(isTRUE(composite)){
    if(isTRUE(stationary)){
      data <- simulateY_comp_s(x, x_test, parameters, seed = seed)
    }else{
      data <- simulateY_comp_ns(x, x_test, parameters, seed = seed)
    }
  }else{
    if(isTRUE(stationary)){
      data <- simulateY_noncomp_s(x, x_test, parameters, seed = seed)
    }else{
      data <- simulateY_noncomp_ns(x, x_test, parameters, seed = seed)
    }
  }

  return(data)
}


simulateY_noncomp_ns <- function(x, x_test, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  n_test <- nrow(x_test)

  R <- get_cor_mat_R(rbind(x, x_test), parameters$rho)
  K <- get_cov_mat_s_R(parameters$sigma2V,
                       R = get_cor_mat_R(rbind(x, x_test), parameters$rhoV),
                       1e-10)
  VAndVTest <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + n_test), K))
  C <- get_cov_mat_ns_R(VAndVTest, R, parameters$sigma2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sigma2eps

  YAndy_test <- MASS::mvrnorm(1, rep(parameters$beta0, n + n_test), C)
  parameters$V <- VAndVTest[1:n]
  parameters$VTest <- VAndVTest[-(1:n)]

  data <- list(y = YAndy_test[1:n],
               y_test = YAndy_test[-(1:n)],
               parameters = parameters)

  return(data)
}

simulateY_noncomp_s <- function(x, x_test, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  n_test <- nrow(x_test)

  R <- get_cor_mat_R(rbind(x, x_test), parameters$rho)
  C <- get_cov_mat_s_R(parameters$sigma2, R, parameters$sigma2eps)
  diag(C)[-(1:n)] <- diag(C)[-(1:n)] - parameters$sigma2eps

  YAndy_test <- MASS::mvrnorm(1, rep(parameters$beta0, n + n_test), C)

  data <- list(y = YAndy_test[1:n],
               y_test = YAndy_test[-(1:n)],
               parameters = parameters)

  return(data)

}


simulateY_comp_ns <- function(x, x_test, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  n_test <- nrow(x_test)

  G <- get_cor_mat_R(rbind(x, x_test), parameters$rhoG)
  L <- get_cor_mat_R(rbind(x, x_test), parameters$rhoL)

  K <- get_cov_mat_s_R(parameters$sigma2V,
                       R = get_cor_mat_R(rbind(x, x_test), parameters$rhoV),
                       1e-10)
  VAndVTest <- exp(MASS::mvrnorm(1, parameters$muV*rep(1, n + n_test), K))

  CG <- parameters$w*get_cov_mat_ns_R(VAndVTest, G, 0)
  CL <- (1 - parameters$w)*get_cov_mat_ns_R(VAndVTest, L, 0)
  CE <- diag(c(rep(parameters$sigma2eps, n), rep(0, n_test)))

  GAndGTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + n_test), CG)
  LAndLTest <- MASS::mvrnorm(1, rep(0, n + n_test), CL)
  EAndETest <- MASS::mvrnorm(1, rep(0, n + n_test), CE)
  YAndy_test <- GAndGTest + LAndLTest + EAndETest

  parameters$V <- VAndVTest[1:n]
  parameters$VTest <- VAndVTest[-(1:n)]

  data <- list(y = YAndy_test[1:n],
               y_test = YAndy_test[-(1:n)],
               yG = GAndGTest[1:n],
               yGTest = GAndGTest[-(1:n)],
               yL = LAndLTest[1:n],
               yLTest = LAndLTest[-(1:n)],
               yE = EAndETest[1:n],
               yETest = EAndETest[-(1:n)],
               parameters = parameters)

  return(data)

}

simulateY_comp_s <- function(x, x_test, parameters, seed){

  set.seed(seed)

  n <- nrow(x)
  n_test <- nrow(x_test)

  G <- get_cor_mat_R(rbind(x, x_test), parameters$rhoG)
  L <- get_cor_mat_R(rbind(x, x_test), parameters$rhoL)

  CG <- parameters$w*get_cov_mat_s_R(parameters$sigma2, G, 0)
  CL <- (1 - parameters$w)*get_cov_mat_s_R(parameters$sigma2, L, 0)
  CE <- diag(c(rep(parameters$sigma2eps, n), rep(0, n_test)))

  GAndGTest <- MASS::mvrnorm(1, rep(parameters$beta0, n + n_test), CG)
  LAndLTest <- MASS::mvrnorm(1, rep(0, n + n_test), CL)
  EAndETest <- MASS::mvrnorm(1, rep(0, n + n_test), CE)
  YAndy_test <- GAndGTest + LAndLTest + EAndETest

  data <- list(y = YAndy_test[1:n],
               y_test = YAndy_test[-(1:n)],
               yG = GAndGTest[1:n],
               yGTest = GAndGTest[-(1:n)],
               yL = LAndLTest[1:n],
               yLTest = LAndLTest[-(1:n)],
               yE = EAndETest[1:n],
               yETest = EAndETest[-(1:n)],
               parameters = parameters)

  return(data)

}

sims_unpack_dots <- function(d, ...){

  dots <- list(...)

  dots_names <- names(dots)
  grid_test <- ifelse("grid_test" %in% dots_names &&
                        (isTRUE(dots$grid_test) || isFALSE(dots$grid_test)),
                      dots$grid_test, FALSE)

  if(grid_test){
    if("grid_test_size" %in% dots_names && is.numeric(dots$grid_test_size) &&
       dots$grid_test_size >= 1){
      grid_test_size <- as.integer(dots$grid_test_size)
    }else{

      if(d <= 4){
        grid_test_size <- switch(d,
                                 100L,
                                 10L,
                                 7L,
                                 5L)

      }else{
        message(strwrap(prefix = " ", initial = "",
                        "A grid in high dimensions quickly gets very large.
                        Please input a smaller 'grid_test_size' or set
                        'grid_test' to FALSE. Proceeding with grid_test =
                        FALSE."))
        grid_test <- FALSE
        grid_test_size <- numeric()
      }
    }
  }else{
    grid_test_size <- numeric()
  }


  randomX1D <- ifelse("randomX1D" %in% dots_names &&
                      (isTRUE(dots$randomX1D) || isFALSE(dots$randomX1D)),
                    dots$randomX1D, FALSE)

  return(list(randomX1D = randomX1D, grid_test = grid_test,
              grid_test_size = grid_test_size))

}
