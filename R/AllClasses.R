## Contains all classes for the bcgp package

#' An S4 class to represent a simulated model
#'
#' This class contains training and test data simulated directly from the BCGP
#' model along with some other information about the type of model the data was
#' simulated from (composite/non-composite, stationary/non-stationary,
#' deterministic/noisy). It will most commonly be used after a call to
#' \code{simulate_from_model()} or \code{bcgpsims()}. The \code{bcgpsims} model
#' can then be plotted or sent to the function \code{bcgp_sampling} to fit the
#' model.
#'
#' @slot training A list with elements \code{x}, an \code{n x d} matrix
#' containing the simulated training data locations, and \code{y}, a vector of
#' length \code{n} containing the simulated training data values. If the
#' \code{composite} argument was \code{TRUE} and the \code{decomposition}
#' argument was \code{TRUE} in the call to \code{simulate_from_model},
#' \code{training} will also have elements \code{yG},
#' \code{yL}, and \code{yE} that correspond to the \strong{G}lobal,
#' \strong{L}ocal, and \strong{E}rror process values at the training data
#' locations.
#' @slot test A list with elements \code{x}, an \code{n_test x d} matrix
#' containing the simulated test data locations, and \code{y}, a vector of
#' length \code{n_test} containing the simulated test data values. If
#' \code{composite} was \code{TRUE} and \code{decomposition} was \code{TRUE} in
#' the call to \code{simulate_from_model}, \code{test} will also have elements
#' \code{yG}, \code{yL}, and \code{yE} that correspond to the \strong{G}lobal,
#' \strong{L}ocal, and \strong{E}rror process values at the test data locations.
#' Note: \code{yE} in \code{test} will \emph{always} be a zero vector. It is
#' only returned for completeness.
#' @slot parameters A list containing the parameters used to simulate from the
#' model.
#' @slot stationary A logical indicating whether the model is stationary or not.
#' @slot composite A logical indicating whether the model is composite or not.
#' @slot seed The seed used by the random number generator to generate the data.
#' @seealso \code{\link{simulate_from_model}} \code{\link{bcgpsims}}
#' @examples
#' simulate_from_model(composite = TRUE, stationary = FALSE, noise = FALSE)
#' bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' @export
setClass(Class = "bcgpsims",
         slots = c(training = "list",
                   test = "list",
                   parameters = "list",
                   stationary = "logical",
                   composite = "logical",
                   seed = "integer"),
         prototype = list(training = list(),
                          test = list(),
                          parameters = list(),
                          stationary = logical(),
                          composite = logical(),
                          seed = NA_integer_))

#' An S4 class to represent a BCGP model
#'
#' This class contains the data, information about the distributions for the
#' parameters, the values of the hyperparameters, along with some other
#' information about the type of model desired to be fit
#' (composite/non-composite, stationary/non-stationary, deterministic/noisy).
#'
#' @slot data A list that contains the training data. One element of the list
#' contains the raw data, and one element contains the scaled data in which the
#' independent variables are scaled to \eqn{[0, 1]^d}, and the response variable
#' is scaled to have mean 0 and variance 1.
#' @slot composite A logical indicating whether the model is composite or not.
#' @slot stationary A logical indicating whether the model is stationary or not.
#' @slot noise A logical indicating whether the data is noisy or deterministic
#' (as from a computer model).
#' @slot priors A list with an element for each parameter in the BCGP model
#' specified by \code{composite} and \code{stationary}. Each element contains
#' the values of the hyperparameters for each parameter.
#' @slot distributions A list with an element for each parameter in the BCGP
#' model specified by \code{composite} and \code{stationary}. Each element
#' contains a character string identifying the prior distribution for each
#' parameter.
#' @seealso \code{\link{bcgpmodel}}
#' @examples
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' bcgpmodel(x = simData@training$x, y = simData@training$y
#'           composite = TRUE, stationary = FALSE, noise = TRUE)
#' @export
setClass(Class = "bcgpmodel",
         slots = c(data = "list",
                   composite = "logical",
                   stationary = "logical",
                   noise = "logical",
                   priors = "list",
                   distributions = "list"),
         prototype = list(data = list(),
                          composite = logical(),
                          stationary = logical(),
                          noise = logical(),
                          priors = list(),
                          distributions = list()))


setOldClass("mcmc.list")

#' An S4 class to represent a fitted BCGP model
#'
#' This class contains the data, both raw and scaled, information about the
#' distributions for the parameters, the values of the hyperparameters, the
#' number of chains desired to be fit, whether the data was scaled to
#' \eqn{[0, 1]^d} before fitting,
#' along with some other information about the type of model fit
#' (composite/non-composite, stationary/non-stationary, deterministic/noisy).
#'
#' @slot model_name A character describing the type of model that was fit.
#' Corresponds to the values for \code{stationary} and \code{composite}.
#' @slot data A list that contains the training data. One element of the list
#' contains the raw data, and one element contains the scaled data in which the
#' independent variables are scaled to \eqn{[0, 1]^d}, and the response variable
#' is scaled to have mean 0 and variance 1.
#' @slot composite A logical indicating whether the model is composite or not.
#' @slot stationary A logical indicating whether the model is stationary or not.
#' @slot noise A logical indicating whether the data is noisy or deterministic
#' (as from a computer model).
#' @slot scaled A logical indicating whether the data was scaled before fitting.
#'  It is highly recommended to scale the data before fitting.
#' @slot chains A positive integer specifying the number of Markov chains
#' @slot priors A list with an element for each parameter in the BCGP model
#' specified by \code{composite} and \code{stationary}. Each element contains
#' the values of the hyperparameters for each parameter.
#' @slot distributions A list with an element for each parameter in the BCGP
#' model specified by \code{composite} and \code{stationary}. Each element
#' contains a character string identifying the prior distribution for each
#' parameter.
#' @slot init A list of length \code{chains} that contains initial values for
#' each parameter for the MCMC algorithm.
#' @slot model_pars A character vector giving the parameter names
#' @slot par_dims A list giving the dimension of each parameter
#' @slot sim A list containing the posterior samples, along with other
#' information
#' @slot algorithm Either \code{"NUTS"} for the \emph{No U-Turn Sampler}
#' implemented by Stan, or \code{"MH"} for a Metropolis-Hastings algorithm
#' @slot sampler_args A list containing information about the sampling
#' algorithm.
#' @seealso \code{\link{bcgp_sampling}} \linkS4class{bcgpmodel}
#' \code{\link{bcgpmodel}}
#' @examples
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' model <- bcgpmodel(x = simData@training$x, y = simData@training$y
#'                    composite = TRUE, stationary = FALSE, noise = TRUE)
#' \dontrun{
#' bcgp_sampling(model, algorithm = "NUTS", scaled = TRUE, chains = 3L,
#'               cores = 1L, control = list(adapt_delta = 0.99,
#'                                          max_treedepth = 15))
#' }
#' @export
setClass(Class = "bcgpfit",
         slots = c(model_name = "character",
                   data = "list",
                   composite = "logical",
                   stationary = "logical",
                   noise = "logical",
                   scaled = "logical",
                   chains = "integer",
                   priors = "list",
                   distributions = "list",
                   init = "list",
                   model_pars = "character",
                   par_dims = "list",
                   sims = "mcmc.list",
                   algorithm = "character",
                   sampler_args = "list",
                   date = "character"))

#' An S4 class to represent a fitted BCGP model with predictions
#'
#' This class contains the data, both raw and scaled, information about the
#' distributions for the parameters, the values of the hyperparameters, the
#' number of chains desired to be fit, whether the data was scaled to
#' \eqn{[0, 1]^d} before fitting,
#' along with some other information about the type of model fit
#' (composite/non-composite, stationary/non-stationary, deterministic/noisy),
#' and predictions at new data locations.
#'
#' @slot model_name A character describing the type of model that was fit.
#' Corresponds to the values for \code{stationary} and \code{composite}.
#' @slot data A list that contains the training data. One element of the list
#' contains the raw data, and one element contains the scaled data in which the
#' independent variables are scaled to \eqn{[0, 1]^d}, and the response variable
#' is scaled to have mean 0 and variance 1.
#' @slot composite A logical indicating whether the model is composite or not.
#' @slot stationary A logical indicating whether the model is stationary or not.
#' @slot noise A logical indicating whether the data is noisy or deterministic
#' (as from a computer model).
#' @slot scaled A logical indicating whether the data was scaled before fitting.
#'  It is highly recommended to scale the data before fitting.
#' @slot chains A positive integer specifying the number of Markov chains
#' @slot priors A list with an element for each parameter in the BCGP model
#' specified by \code{composite} and \code{stationary}. Each element contains
#' the values of the hyperparameters for each parameter.
#' @slot distributions A list with an element for each parameter in the BCGP
#' model specified by \code{composite} and \code{stationary}. Each element
#' contains a character string identifying the prior distribution for each
#' parameter.
#' @slot init A list of length \code{chains} that contains initial values for
#' each parameter for the MCMC algorithm.
#' @slot model_pars A character vector giving the parameter names
#' @slot par_dims A list giving the dimension of each parameter
#' @slot sim A list containing the posterior samples, along with other
#' information
#' @slot algorithm Either \code{"NUTS"} for the \emph{No U-Turn Sampler}
#' implemented by Stan, or \code{"MH"} for a Metropolis-Hastings algorithm
#' @slot sampler_args A list containing information about the sampling
#' algorithm.
#' @slot preds A list containing the prediction locations and the predictions.
#' @seealso \code{\link{bcgp_sampling}} \linkS4class{bcgpmodel}
#' \linkS4class{bcgpfit} \code{\link[bcgp]{predict}}
#' \code{\link[bcgp]{posterior_predict}}
#' @examples
#' simData <- bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#' model <- bcgpmodel(x = simData@training$x, y = simData@training$y
#'                    composite = TRUE, stationary = FALSE, noise = TRUE)
#' \dontrun{
#' bcgp_sampling(model, algorithm = "NUTS", scaled = TRUE, chains = 3L,
#'               cores = 1L, control = list(adapt_delta = 0.99,
#'                                          max_treedepth = 15))
#' }
#' @export
setClass(Class = "bcgpfitpred",
         slots = c(preds = "list"),
         contains = "bcgpfit")
