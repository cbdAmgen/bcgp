% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcgpfitpred-class.R
\name{predictive_interval,bcgpfitpred-method}
\alias{predictive_interval,bcgpfitpred-method}
\title{\code{predictive_interval} method for a \code{bcgpfitpred} object}
\usage{
\S4method{predictive_interval}{bcgpfitpred}(object, prob = 0.95, type = c("central", "HPD"))
}
\arguments{
\item{object}{a bcgpfitpred object}

\item{prob}{A numeric scalar between 0 and 1 specifying the uncertainty
interval width. Defaults to 0.95.}

\item{type}{A character specifying whether "central" intervals based on
quantiles or highest posterior density ("HPD") intervals should be
constructed. Defaults to "central".}
}
\value{
A matrix with two columns and as many rows as the number of
prediction locations For a given value of \code{prob}, \eqn{p}, the columns
correspond to the lower and upper \eqn{100p}\% interval limits and have the
names \code{lower} and \code{upper}.
}
\description{
This calculates Bayesian predictive intervals
}
\examples{
data_sim <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE,
                     d = 2, decomposition = TRUE)

model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
                   composite = TRUE, stationary = FALSE, noise = FALSE)
fit <- bcgp_sampling_with_predictions(model, x_pred = data_sim@test$x,
                                      scaled = FALSE, cores = 4,
                                      n_mcmc = 500, burnin = 200)

predictive_interval(fit, prob = 0.90, type = "central")

}
