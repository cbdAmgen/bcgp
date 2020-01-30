#' Create an object of class bcgpsims
#'
#' This function creates an object of class bcgpsims, simulating data from the
#' specified (by stationary, composite, noise) Gaussian process model.
#'
#' \code{bcgpsims} returns an instance of S4 class \code{bcgpsims}.
#' This object can then be plotted to get an idea of what draws from these
#' models look like, or the data in the object can be fit.
#' \code{simulate_from_model()} can also be called to create a \code{bcgpsims}
#' object that contains simulated data.
#'
#' @param composite A logical, \code{TRUE} for a composite of a global process
#' and a local process, \code{FALSE} for non-composite. Defaults to \code{TRUE}.
#' @param stationary A logical, \code{FALSE} for a non-stationary process,
#' \code{TRUE} for a stationary process. If \code{FALSE}, the variance for the
#' process is \eqn{\sigma^2(x)}, and if \code{TRUE}, the variance is
#' \eqn{\sigma^2}. Defaults to \code{FALSE}.
#' @param noise If the data should be noise-free (such as from a deterministic
#' computer model), then \code{noise} should be \code{FALSE}. Otherwise, it
#' should be \code{TRUE}. Defaults to \code{FALSE}.
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
#'   \item{\code{grid_test}}{A logical indicating whether the test data should be
#'   generated on a grid. Defaults to \code{FALSE}. Only useful for
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
#' @seealso \linkS4class{bcgpsims} \code{\link{simulate_from_model}}
#' \code{\link{create_parameter_list}}
#' @examples
#' bcgpsims(composite = TRUE, stationary = FALSE, noise = FALSE)
#'
#' params <- create_parameter_list()
#' params$w <- 0.99
#' bcgpsims(parameters = params, randomX1D = TRUE)
#' @export
bcgpsims <- function(composite = TRUE, stationary = FALSE,
                     noise = FALSE, d = 1L, n = 15*d, n_test = 100*d,
                     parameters = create_parameter_list(composite,
                                                        stationary,
                                                        noise, d),
                     seed = sample.int(.Machine$integer.max, 1),
                     ...) {
  simulate_from_model(composite, stationary, noise,
                      d, n, n_test,
                      parameters,
                      seed,
                      ...)
}

#' Plot a bcgpsims object
#'
#' This plots the simulated data in a \code{bcgpsims} object.
#'
#' Plotting simulated data helps to understand what these types of models look
#' like and understand the interaction between the variance process
#' \eqn{\sigma^2(x)} and the data process \eqn{Y(x)} for non-stationary models.
#' Plotting is only supported for one and two-dimensional data.
#'
#' @param x An instance of class \linkS4class{bcgpsims}
#' @param ... optional parameters
#' \describe{
#'   \item{\code{raster}}{A logical indicating whether the 2-dimensional process
#'   should be plotted in points or smoothed out as in
#'   \code{\link[ggplot2]{geom_raster}}. If the \code{bcgpsims} object was
#'   created with randomly generated test data (i.e., not in a grid), this
#'   argument is ignored, and the plot will be points and not smoothed out.
#'   Only useful for 2-D data.}
#' }
#' @param process Which process to plot: the data process and/or the variance
#' process
#' @param decomposition A logical indicating whether to plot the data process
#' decomposed into global and local processes (TRUE) or the overall process by
#' itself (FALSE). This parameter is ignored when the data has dimension larger
#' than 1. Defaults to FALSE
#' @param print A logical indicating whether to automatically print the plots.
#' Defaults to TRUE.
#' @return Either one \code{\link[ggplot2]{ggplot}} object or a list of ggplot
#' objects that can be further customized using the \pkg{ggplot2} package.
#' @seealso \linkS4class{bcgpsims}
#' @examples
#' nsComp <- simulate_from_model(composite = TRUE, stationary = FALSE, d = 1,
#'                               decomposition = TRUE)
#' z <- plot(nsComp, print = FALSE, decomposition = TRUE)
#'
#' sNonComp <- simulate_from_model(composite = FALSE, stationary = TRUE, d = 2,
#'                                 grid_test = TRUE, grid_test_size = 10)
#' plot(sNonComp, process = "y", raster = TRUE)
#' plot(sNonComp, process = "y", raster = FALSE)
#' @export
setMethod("plot", signature(x = "bcgpsims"),
          function(x, ..., process = c("y", "variance"),
                   decomposition = FALSE, print = TRUE){
            d <- ncol(x@training$x)
            process <- match.arg(process, several.ok = TRUE)

            if("y" %in% process){
              yPlot <- plotDataSims(x, decomposition, ...)
              if(print) print(yPlot)
              if(length(process) == 1) return(yPlot)
            }

            if("variance" %in% process){
              if(x@stationary){
                warning(strwrap(prefix = " ", initial = "",
                                "No need to plot the variance process, since the
                                variance is constant in a stationary process."))
                varPlot <- ggplot2::ggplot()
              }else{
                varPlot <- plotVarSims(x, ...)
                if(print) print(varPlot)
                if(length(process) == 1) return(varPlot)
              }
            }
            return(invisible(list(dataProcess = yPlot,
                                  varianceProcess = varPlot)))
          })
