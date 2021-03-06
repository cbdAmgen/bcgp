% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcgpmodel-class.R
\name{bcgp_sampling_and_prediction,bcgpmodel-method}
\alias{bcgp_sampling_and_prediction,bcgpmodel-method}
\title{Draw samples from a BCGP model and make predictions}
\usage{
\S4method{bcgp_sampling_and_prediction}{bcgpmodel}(
  object,
  x_pred,
  algorithm = c("NUTS", "MH"),
  scaled = TRUE,
  chains = 4L,
  cores = getOption("mc.cores", 1L),
  init = "random",
  num_updates = 5,
  num_adapt = 1000,
  burnin = 1000,
  n_mcmc = 10000,
  thin = 1,
  ...,
  control = NULL
)
}
\arguments{
\item{object}{An instance of S4 class \code{bcgpmodel}.}

\item{x_pred}{A matrix with the same number of columns as the data matrix in
\code{object}.}

\item{algorithm}{Either \code{"NUTS"} for the \code{No U-Turn Sampler}
implemented by Stan, or \code{"MH"} for a Metropolis-Hastings algorithm}

\item{scaled}{A logical indicating whether the data should be scaled before
fitting. It is highly recommended to scale the data before fitting. Defaults
to \code{TRUE}. Currently, only \code{scaled = TRUE} is supported, but that
will change in the near future.}

\item{chains}{A positive integer specifying the number of Markov chains}

\item{cores}{The number of cores to use when executing the Markov chains in
parallel. The default is to use the value of the \code{mc.cores} option if it
has been set and otherwise to default to 1 core. However, it is recommended
to set it to be as many processors as the hardware and RAM allow (up to the
number of chains). See \code{\link[parallel]{detectCores}}if you don't know
this number for your system.}

\item{init}{Can be either the string "random" or a list of length
\code{chains}. The elements of this list should be named lists, where each of
these has the name of a parameter and is used to specify the initial values
for that parameter for the corresponding chain.

\describe{
\code{init = "random"} (default):
The initial values will be generated randomly from their respective prior
distributions.

\code{init} via list:
Set initial values by providing a list equal in length to the number of
Markov chains. A call to \code{create_inits()} will assist in the correct
creation of this list. This is not currently supported, but will be in the
future.
}}

\item{num_updates}{A positive integer for the number of updates in the
proposal stepsize adaptation phase. Ignored if \code{algorithm = "NUTS"}.}

\item{num_adapt}{A positive integer for the number of samples within each
update in the proposal stepsize adaptation phase. Ignored if
\code{algorithm = "NUTS"}.}

\item{burnin}{A positive integer for the number of burnin samples to discard
after the stepsize adaptation phase is finished. This is equivalent to the
parameter \code{warmup} in \code{\link[rstan]{stan}}.}

\item{n_mcmc}{The number of samples to be kept for each Markov chain.}

\item{thin}{A positive integer specifying the period for saving samples. The
default is 1, and this number should not be changed, as thinning isn't
necessary in these models, and it throws away information. Currently, only
\code{thin = 1} is supported, and this argument may be deprecated in the
future.}

\item{...}{optional parameters, only if \code{algorithm = "NUTS"}. See the
documentation for \code{\link[rstan]{stan}}.}

\item{control}{A named list of parameters to control the NUTS algorithm's
behavior. It defaults to NULL so all the default values are used. Ignored
unless  \code{algorithm = "NUTS"}.}
}
\value{
An instance of S4 class \code{bcgpfitpred}. It contains information
about the data, model, sampling algorithm, and sample draws from
the posterior for the parameters and new predictions.
}
\description{
\code{bcgp_sampling_and_prediction} draws samples from the model defined by
class \code{bcgpmodel} and predicts at the locations in \code{x_pred}
}
\details{
This returns an instance of S4 class \code{bcgpfitpred}. It contains
information about the data, model, sampling algorithm, and sample draws from
the posterior for the parameters and new predictions.
}
