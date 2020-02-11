setGeneric(name = "print_predictions",
           def = function(object, ...) {
             standardGeneric("print_predictions")
           })

#' @export
setMethod("print_predictions", "bcgpfitpred",
          function(object, pars = c("y", "V"),
                   quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                   digits_summary = 2) {


            pars <- match.arg(pars, several.ok = TRUE)

            summ <- summary_predictions(object, pars, quantiles)

            n_mcmc <- object@sampler_args$general$n_mcmc
            burnin <- object@sampler_args$general$burnin
            thin <- object@sampler_args$general$thin

            cat("Inference for BCGP model: ", object@model_name, '.\n',
                sep = '')
            cat(object@chains, " chains, each with n_mcmc = ", n_mcmc,
                "; burnin = ", burnin, "; thin = ", thin, "; \n",
                "total post-burnin draws = ", object@chains*n_mcmc, ".\n\n",
                sep = '')

            print(round(summ, digits_summary))

            cat("\nSamples were drawn using ", object@algorithm, " at ",
                object@date, ".\n",
                "For each parameter, n_eff is a crude measure of effective sample\n",
                "size, and Rhat is the potential scale reduction factor on split\n",
                "chains (at convergence, Rhat = 1).\n", sep = '')

            return(invisible(NULL))

          })



setGeneric(name = "summary_predictions",
           def = function(object, ...) {
             standardGeneric("summary_predictions")
           })

#' \code{summary_predictions} method for a \code{bcgpfitpred} object
#'
#' This gives a summary of the posterior predictions contained in a
#' \code{bcgpfitpred} object
#'
#' @param object a bcgpfitpred object
#' @param pars a character vector specifying the parameters to summarize. One or
#' both of "y" and "V" may be specified to summarize the posterior predictions
#' of the process or the variance process, respectively. If the process is
#' stationary, then thie argument will be ignored.
#' @param quantiles a numeric vector specifying the desired quantiles
#'
#' @examples
#' data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                      d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling_with_predictions(model, x_pred = data_sim@test$x,
#'                                       scaled = FALSE, cores = 4,
#'                                       n_mcmc = 500, burnin = 200)
#'
#' fit
#' print(fit, pars = c("beta0", "w", "rhoG", "rhoL"), digits_summary = 3)
#' print_predictions(fit, pars = c("y", "V"), digits_summary = 3)
#' summary(fit)
#' summary_predictions(fit)
#'
#' @export
setMethod("summary_predictions", signature = "bcgpfitpred",
          function(object, pars = c("y", "V"),
                   quantiles = c(0.025, 0.25, 0.50, 0.75, 0.975), ...) {

            pars <- match.arg(pars, several.ok = TRUE)

            if("V" %in% pars && isTRUE(object@stationary)){

              warning(strwrap(prefix = " ", initial = "",
                              "The variance process is constant for a
                              stationary model."))
              if("y" %in% pars){
                pars <- "y"
                warning(strwrap(prefix = " ", initial = "",
                                "Only the process predictions will be
                                printed."))
              }else{
                stop(strwrap(prefix = " ", initial = "",
                             "For a stationary process, please only print
                             'y_pred'."))
              }
            }
            out <- do.call(rbind, lapply(object@preds[pars], coda_summary,
                                         quantiles = quantiles))

            return(out)
          }
)

setGeneric(name = "rmspe",
           def = function(object, ...) {
             standardGeneric("rmspe")
           })

#' \code{rmspe} method for a \code{bcgpfitpred} object
#'
#' This calculates the root mean square prediction error for the posterior
#' predictions contained in a \code{bcgpfitpred} object in comparison to a
#' user-provided vector
#'
#' @param object a bcgpfitpred object
#' @param truth a numeric vector specifying the true values that should be
#' compared to the \code{bcgp} predictions
#' @param pred_type a character specifying the posterior point estimator to use.
#' It can be either "mean" or "median". Defaults to "mean
#'
#' @examples
#' data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                      d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling_with_predictions(model, x_pred = data_sim@test$x,
#'                                       scaled = FALSE, cores = 4,
#'                                       n_mcmc = 500, burnin = 200)
#'
#' rmspe(fit, truth = data_sim@test$y, pred_type = "mean")
#'
#' @export
setMethod("rmspe", signature = "bcgpfitpred",
          function(object, truth, pred_type = c("mean", "median")) {

            check_valid_truth(object, truth)
            pred_type <- match.arg(pred_type)

            if(coda::is.mcmc.list(object@preds$y)){
              y_matrix <- as.matrix(object@preds$y)
            }else{
              stop(strwrap(prefix = " ", initial = "",
                           "The fitted object with predictions is corrupted."))

            }
            pred <- apply(y_matrix, 2, pred_type)

            sqrt(mean((truth - pred)^2))

          }
)


setGeneric(name = "plot_predictions",
           def = function(object, ...) {
             standardGeneric("plot_predictions")
           })

#' \code{plot_predictions} method for a \code{bcgpfitpred} object
#'
#' This plots the predictions in a \code{bcgpfitpred} object.
#'
#' Plotting the process \eqn{Y(x)}, helps to visualize the predictions.
#' Plotting is only supported for one-dimensional data.
#'
#' @param object An instance of class \linkS4class{bcgpfitpred}
#' @param prob A numeric scalar between 0 and 1 specifying the uncertainty
#' interval width. Defaults to 0.95. \code{NULL} will give no uncertainty
#' intervals.
#' @param estimate A character specifying whether the mean or median should be
#' used as the point estimate. Defaults to "mean".
#' @param interval_type A character specifying whether "central" intervals based
#' on quantiles or highest posterior density ("HPD") intervals should be
#' constructed. Defaults to "central". If \code{prob} is \code{NULL} then this
#' argument is ignored.
#' @param overlay_data A logical specifying whether the observed data should be
#' overlaid. Defaults to TRUE.
#' @param print A logical indicating whether to automatically print the plots.
#' Defaults to TRUE.
#' @return A \code{\link[ggplot2]{ggplot}} object that can be further customized
#'  using the \pkg{ggplot2} package.
#' @seealso \linkS4class{bcgpfitpred}
#' @examples
#' data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                      d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling_with_predictions(model, x_pred = data_sim@test$x,
#'                                       scaled = FALSE, cores = 4,
#'                                       n_mcmc = 500, burnin = 200)
#'
#' plot(fit, prob = 0.90, esimate = "mean", overlay_data = TRUE,
#'      interval_type = "central")
#' @export
setMethod("plot_predictions", signature(object = "bcgpfitpred"),
          function(object, prob = 0.95, estimate = c("mean", "median"),
                   interval_type = c("central", "HPD"), overlay_data = TRUE,
                   print = TRUE){

            if(ncol(object@data$raw$x) != 1)
              stop(strwrap(prefix = " ", initial = "",
                           "Currently plotting is only supported for
                            one-dimensional data."))
            if(!is.null(prob)) check_prob(prob)

            estimate <- match.arg(estimate)
            interval_type <- match.arg(interval_type)

            if(coda::is.mcmc.list(object@preds$y)){
              predictions_matrix <- as.matrix(object@preds$y)
            }else{
              stop(strwrap(prefix = " ", initial = "",
                           "The fitted object with predictions is corrupted."))

            }

            point_estimates <- matrix(apply(predictions_matrix, 2, estimate),
                                      ncol = 1)
            colnames(point_estimates) <- "estimate"

            if(!is.null(prob)){
              interval_estimates <-
                as.data.frame(predictive_interval(object, prob = prob,
                                                  type = interval_type))
            }else{
              interval_estimates <- data.frame(lower = NA, upper = NA)
            }


            predictions_df <- cbind(data.frame(x = object@preds$x),
                                    point_estimates, interval_estimates)

            data_df <- data.frame(x = object@data$raw$x, y = object@data$raw$y)


            title_stat <- ifelse(object@stationary, "Stationary",
                                 "Non-Stationary")
            title_comp <- ifelse(object@composite, "Composite", "Non-Composite")
            plot_title <- paste("BCGP Predictions for Y(x):", title_stat,
                                title_comp, sep = " ")

            plot_out <- ggplot2::ggplot(data = predictions_df) +
              ggplot2::geom_line(ggplot2::aes(x = x, y = estimate), size = 1.5,
                                 color = "green") +
              ggplot2::theme_classic() +
              ggplot2::ylab("Y(x)") +
              ggplot2::ggtitle("Predictions for Y(x)") +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

            if(!is.null(prob)){
              plot_out <- plot_out +
                ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = lower,
                                                  ymax = upper),
                                     fill = "yellow", alpha = 0.25)
            }

            if(isTRUE(overlay_data)){
              plot_out <- plot_out +
                ggplot2::geom_point(data = data_df,
                                    mapping = ggplot2::aes(x = x, y = y),
                                    color = "red")
            }

            if(print) print(plot_out)
            return(invisible(plot_out))
          })


setGeneric(name = "predictive_interval",
           def = function(object, ...) {
             standardGeneric("predictive_interval")
           })

#' \code{predictive_interval} method for a \code{bcgpfitpred} object
#'
#' This calculates Bayesian predictive intervals
#'
#' @param object a bcgpfitpred object
#' @param prob A numeric scalar between 0 and 1 specifying the uncertainty
#' interval width. Defaults to 0.95.
#' @param type A character specifying whether "central" intervals based on
#' quantiles or highest posterior density ("HPD") intervals should be
#' constructed. Defaults to "central".
#' @return A matrix with two columns and as many rows as the number of
#' prediction locations.
#'
#' @examples
#' data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
#'                      d = 2, decomposition = TRUE)
#'
#' model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
#'                    composite = TRUE, stationary = FALSE, noise = FALSE)
#' fit <- bcgp_sampling_with_predictions(model, x_pred = data_sim@test$x,
#'                                       scaled = FALSE, cores = 4,
#'                                       n_mcmc = 500, burnin = 200)
#'
#' predictive_interval(fit, prob = 0.90, type = "central")
#'
#' @export
setMethod("predictive_interval", signature = "bcgpfitpred",
          function(object, prob = 0.95, type = c("central", "HPD")) {

            check_prob(prob)
            type <- match.arg(type)

            if(coda::is.mcmc.list(object@preds$y)){

              prediction_matrix <- as.matrix(object@preds$y)

              if(type == "central"){
                out <- central_intervals(prediction_matrix, prob)
              }else{
                out <- hpd_intervals(prediction_matrix, prob)
              }
            }else{
              stop(strwrap(prefix = " ", initial = "",
                           "The fitted object with predictions is corrupted."))

            }

            return(out)
          })
