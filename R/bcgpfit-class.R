#' @export
setMethod("show", signature = "bcgpfit",
          function(object){
            print.bcgpfit(object, pars = object@model_pars)

          })

##' @method print bcgpfit
##' @export
print.bcgpfit <- function(x, pars = x@model_pars,
                          quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                          digits_summary = 2){

  summ <- summary(x, pars, quantiles)

  n_mcmc <- x@sampler_args$general$n_mcmc
  burnin <- x@sampler_args$general$burnin
  thin <- x@sampler_args$general$thin

  cat("Inference for BCGP model: ", x@model_name, '.\n', sep = '')
  cat(x@chains, " chains, each with n_mcmc = ", n_mcmc,
      "; burnin = ", burnin, "; thin = ", thin, "; \n",
      "total post-burnin draws = ", x@chains*n_mcmc, ".\n\n", sep = '')

  print(round(summ, digits_summary))

  cat("\nSamples were drawn using ", x@algorithm, " at ", x@date, ".\n",
      "For each parameter, n_eff is a crude measure of effective sample\n",
      "size, and Rhat is the potential scale reduction factor on split\n",
      "chains (at convergence, Rhat = 1).\n", sep = '')

  return(invisible(NULL))

}

#' \code{summary} method for a \code{bcgpfit} object
#'
#' This gives a summary of the posterior draws contained in a \code{bcgpfit}
#' object
#'
#' @param object a bcgpfit object
#' @param pars a character vector specifying the parameters to summarize
#' @param quantiles a numeric vector specifying the desired quantiles for each
#' parameter
#'
#' @examples
#'
#'
#' data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                      d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling(model, scaled = TRUE, cores = 4, n_mcmc = 500,
#'                      burnin = 200)
#'
#' fit
#' print(fit, pars = c("beta0", "w", "rhoG", "rhoL"), digits_summary = 3)
#' summary(fit)
#'
#' @export
setMethod("summary", signature = "bcgpfit",
          function(object, pars,
                   quantiles = c(0.025, 0.25, 0.50, 0.75, 0.975), ...) {

            samples <- object@sims

            if(missing(pars)) pars <- object@model_pars
            if(missing(quantiles))
              quantiles <- c(0.025, 0.25, 0.50, 0.75, 0.975)

            n_eff <- as.matrix(floor(coda::effectiveSize(samples)))
            colnames(n_eff) <- "n_eff"
            R_hats <- as.matrix(coda::gelman.diag(samples)$psrf[, "Point est."])
            colnames(R_hats) <- "Rhat"
            long_sum <- summary(samples, quantiles)

            allPars <- cbind(long_sum$statistics, long_sum$quantiles,
                             n_eff, R_hats)
            pars2 <- paste0("^", pars)
            out <- allPars[grepl(paste(pars2, collapse = "|"),
                                 row.names(allPars)),]

            return(out)

          }
)

setGeneric(name = "posterior_interval",
           def = function(object, ...) {
             standardGeneric("posterior_interval")
           })

#' \code{posterior_interval} method for a \code{bcgpfit} object
#'
#' This calculates Bayesian uncertainty intervals, commonly called
#' \emph{credible} intervals
#'
#' @param object a bcgpfitpred object
#' @param prob A numeric scalar between 0 and 1 specifying the uncertainty
#' interval width. Defaults to 0.95.
#' @param type A character specifying whether "central" intervals based on
#' quantiles or highest posterior density ("HPD") intervals should be
#' constructed. Defaults to "central".
#' @param pars A character vector specifying the parameters for which to return
#' intervals.
#' @return A matrix with two columns and as many rows as model parameters
#' (or a subset of parameters specified by the user). For a given value of
#' \code{prob}, \eqn{p}, the columns correspond to the lower and upper
#' \eqn{100p}\% interval limits and have the names \code{lower} and
#' \code{upper}.
#'
#' @examples
#' data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                      d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling(model, scaled = FALSE, cores = 4,
#'                      n_mcmc = 500, burnin = 200)
#'
#' posterior_interval(fit, prob = 0.90, type = "central",
#'                    pars = c("beta0", "w", "rhoG", "rhoL"))
#'
#' @export
setMethod("posterior_interval", signature = "bcgpfit",
          function(object, prob = 0.95, type = c("central", "HPD"), pars) {


            if(missing(pars)) pars <- object@model_pars

            check_prob(prob)
            type <- match.arg(type)

            if(coda::is.mcmc.list(object@sims)){

              pars2 <- paste0("^", pars)
              pars_to_keep <- grepl(paste(pars2, collapse = "|"),
                                    colnames(object@sims[[1]]))
              parameter_matrix <- as.matrix(object@sims[, pars_to_keep,
                                                        drop = TRUE])

              if(type == "central"){
                out <- central_intervals(parameter_matrix, prob)
              }else{
                out <- hpd_intervals(parameter_matrix, prob)
              }
            }else{
              stop(strwrap(prefix = " ", initial = "",
                           "The fitted object is corrupted."))

            }

            return(out)
          })
