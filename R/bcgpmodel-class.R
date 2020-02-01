#' Create an object of class bcgpmodel
#'
#' \code{bcgpmodel} returns an instance of S4 class \code{bcgpmodel}
#'
#' This creates an instance of S4 class \code{bcgpmodel} that contains
#' user-provided data, user-provided information about the type of model desired
#' to be fit (composite/non-composite, stationary/non-stationary,
#' deterministic/noisy), default values for the hyperparameters for the prior
#' distributions, and information about the prior distributions.
#'
#' @param x An \code{n x d} numeric matrix, where \code{d} is the dimension of
#' the data, and \code{n} is the number of training data observations.
#' @param y A numeric vector of length \code{n}.
#' @param composite A logical, \code{TRUE} for a composite of a global process,
#' a local process, and an error process, \code{FALSE} for non-composite.
#' Defaults to \code{TRUE}.
#' @param stationary A logical, \code{FALSE} for a non-stationary process,
#' \code{TRUE} for a stationary process. If \code{FALSE}, the variance for the
#' process is \eqn{\sigma^2(x)}, and if \code{TRUE}, the variance is
#' \eqn{\sigma^2}. Defaults to \code{FALSE}.
#' @param noise If the data is known to be noise-free (such as from a
#' deterministic computer model), then \code{noise} should be \code{FALSE}.
#' Otherwise, it should be \code{TRUE}. Defaults to \code{FALSE}
#'
#' @return An instance of S4 class \code{bcgpmodel} containing the data, default
#' values for all the prior parameters, information about the process, and
#' information about the distributions.
#' @family preprocessing functions
#' @seealso \linkS4class{bcgpmodel} \code{\link{bcgpfit}}
#' \code{\link{bcgp}}
#' @examples
#' data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' bcgpmodel(x = data_sim@training$x, y = simData@training$y,
#'           composite = TRUE, stationary = FALSE, noise = FALSE)
#' @export

bcgpmodel <- function(x, y, composite = TRUE, stationary = FALSE,
                      noise = FALSE){

  if(is.vector(x)) x <- as.matrix(x)
  validate_bcgpmodel_inputs(x, y, composite, stationary, noise)

  data <- list(x = x, y = y)

  prior_info <- create_priors(composite, stationary, noise, d = ncol(x))

  priors <- prior_info$priors
  distributions <- prior_info$distributions

  # if(scaled){
  #   xIn <- data$scaled$x
  # }else{
  #   xIn <- data$raw$x
  # }
  #
  # inits <- create_inits(x = xIn, composite, stationary, noise,
  #                       priors, chains)


  new("bcgpmodel",
      data = data,
      composite = composite,
      stationary = stationary,
      noise = noise,
      priors = priors,
      distributions = distributions)

}


#' setGeneric(name = "update_inits",
#'            def = function(object, ...) standardGeneric("update_inits"))
#'
#' #' Update the starting values for MCMC chains
#' #'
#' #' \code{update_inits} returns an instance of S4 class \code{bcgpmodel} after
#' #' updating the starting values for each Markov chain
#' #'
#' #' This returns an instance of S4 class \code{bcgpmodel}. It is generally meant
#' #' to update the starting values for the Markov chains after the user has
#' #' changed the prior distributions in an instance of class \code{bcgpmodel}. The
#' #' new initial values are simulated from the prior distribution, ensuring they
#' #' are within the support of the posterior distribution.
#' #'
#' #' @param object An instance of S4 class \code{bcgpmodel}.
#' #' @return An instance of S4 class \code{bcgpmodel} containing initial values
#' #' randomly from the prior distributions for all the parameters
#' #' @family preprocessing functions
#' #' @seealso \linkS4class{bcgpmodel}
#' #' @examples
#' #' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' #' model <- bcgpmodel(x = simData@training$x, y = simData@training$y,
#' #'                    composite = TRUE, stationary = FALSE, noise = FALSE,
#' #'                    scaled = TRUE, chains = 4L)
#' #' model@priors$w$lower <- 0.75
#' #' model <- update_inits(model)
#' #' @export
#'
#' setMethod("update_inits", signature = "bcgpmodel",
#'           function(object){
#'
#'             data <- object@data
#'
#'             if(object@scaled){
#'               xIn <- data$scaled$x
#'             }else{
#'               xIn <- data$raw$x
#'             }
#'
#'             object@inits <- create_inits(x = xIn, object@composite,
#'                                          object@stationary, object@noise,
#'                                          object@priors, object@chains)
#'
#'
#'             return(object)
#'
#'           })



setGeneric(name = "bcgp_sampling",
           def = function(object, ...) { standardGeneric("bcgp_sampling")})


#' Draw samples from a BCGP model
#'
#' \code{bcgp_sampling} draws samples from the model defined by class
#' \code{bcgpmodel}
#'
#' This returns an instance of S4 class \code{bcgpfit}. It contains information
#' about the data, model, sampling algorithm, and sample draws from
#' the posterior
#'
#' @param object An instance of S4 class \code{bcgpmodel}.
#' @param algorithm Either \code{"NUTS"} for the \code{No U-Turn Sampler}
#' implemented by Stan, or \code{"MH"} for a Metropolis-Hastings algorithm
#' @param scaled A logical indicating whether the data should be scaled before
#' fitting. It is highly recommended to scale the data before fitting. Defaults
#' to \code{TRUE}. Currently, only \code{scaled = TRUE} is supported, but that
#' will change in the near future.
#' @param chains A positive integer specifying the number of Markov chains
#' @param cores The number of cores to use when executing the Markov chains in
#' parallel. The default is to use the value of the \code{mc.cores} option if it
#' has been set and otherwise to default to 1 core. However, it is recommended
#' to set it to be as many processors as the hardware and RAM allow (up to the
#' number of chains). See \code{\link[parallel]{detectCores}}if you don't know
#' this number for your system.
#' @param init Can be either the string "random" or a list of length
#' \code{chains}. The elements of this list should be named lists, where each of
#' these has the name of a parameter and is used to specify the initial values
#' for that parameter for the corresponding chain.
#'
#' \describe{
#' \code{init = "random"} (default):
#' The initial values will be generated randomly from their respective prior
#' distributions.
#'
#' \code{init} via list:
#' Set initial values by providing a list equal in length to the number of
#' Markov chains. A call to \code{create_inits()} will assist in the correct
#' creation of this list. This is not currently supported, but will be in the
#' future.
#' }
#' @param num_updates A positive integer for the number of updates in the
#' proposal stepsize adaptation phase. Ignored if \code{algorithm = "NUTS"}.
#' @param num_adapt A positive integer for the number of samples within each
#' update in the proposal stepsize adaptation phase. Ignored if
#' \code{algorithm = "NUTS"}.
#' @param burnin A positive integer for the number of burnin samples to discard
#' after the stepsize adaptation phase is finished. This is equivalent to the
#' parameter \code{warmup} in \code{\link[rstan]{stan}}.
#' @param n_mcmc The number of samples to be kept for each Markov chain.
#' @param thin A positive integer specifying the period for saving samples. The
#' default is 1, and this number should not be changed, as thinning isn't
#' necessary in these models, and it throws away information. Currently, only
#' \code{thin = 1} is supported, and this argument may be deprecated in the
#' future.
#' @param ... optional parameters, only if \code{algorithm = "NUTS"}. See the
#' documentation for \code{\link[rstan]{stan}}.
#' @param control A named list of parameters to control the NUTS algorithm's
#' behavior. It defaults to NULL so all the default values are used. Ignored
#' unless  \code{algorithm = "NUTS"}.
#' @return An instance of S4 class \code{bcgpfit}. It contains information
#' about the data, model, sampling algorithm, and sample draws from
#' the posterior.
#'
#' @export
setMethod("bcgp_sampling", "bcgpmodel",
          function(object, algorithm = c("NUTS", "MH"), scaled = TRUE,
                   chains = 4L, cores = getOption("mc.cores", 1L),
                   init = "random", num_updates = 5, num_adapt = 1000,
                   burnin = 1000, n_mcmc = 10000, thin = 1, ...,
                   control = NULL) {

            ## TODO: error checking for all the inputs here

            if(!isTRUE(scaled)){
              warning(strwrap(prefix = " ", initial = "",
                              "It is currently required that the data be scaled before
                    fitting. Changing 'scaled' to be TRUE."))
              scaled <- TRUE
            }

            if(thin != 1){
              warning(strwrap(prefix = " ", initial = "",
                              "Thinning is not currently supported and is
                              generally unnecessary. Changing 'thin' to 1."))
            }
            thin <- 1

            algorithm <- match.arg(algorithm)

            if(algorithm == "NUTS"){

              warmup <- burnin
              iter <- n_mcmc + warmup
              out1 <- bcgp_stan(object, scaled, chains, cores, iter, warmup,
                                thin, control, ...)

              sampler_args <- list(Stan = out1$sampler_args,
                                   general = list(chains = chains,
                                                  burnin = burnin,
                                                  n_mcmc = n_mcmc,
                                                  thin = thin))

            }
            else{
              warning(strwrap(prefix = " ", initial = "",
                              "Metropolis_Hastings and Gibbs sampler currently
                              not supported. Model will be fit with Stan using
                              the NUTS algorithm."))

              warmup <- burnin
              iter <- n_mcmc + warmup
              out1 <- bcgp_stan(object, scaled, chains, cores, iter, warmup,
                                thin, control, ...)

              sampler_args <- list(Stan = out1$sampler_args,
                                   general = list(chains = chains,
                                                  burnin = burnin,
                                                  n_mcmc = n_mcmc,
                                                  thin = thin))
              # out1 <- bcgp_MH(object, scaled, ...)
              #
              # sampler_args <- list(MH = list(num_updates = 5,
              #                                num_adapt = 1000),
              #                      general = list(chains = chains,
              #                                     burnin = burnin,
              #                                     n_mcmc = n_mcmc,
              #                                     thin = thin))
            }

            # return(out)
            new("bcgpfit",
                model_name = out1$model_name,
                data = out1$data,
                stationary = object@stationary,
                composite = object@composite,
                noise = object@noise,
                scaled = scaled,
                chains = chains,
                priors = object@priors,
                distributions = object@distributions,
                init = out1$init,
                model_pars = out1$model_pars,
                par_dims = out1$par_dims,
                sims = out1$sims,
                algorithm = algorithm,
                sampler_args = sampler_args,
                date = date())

          })
