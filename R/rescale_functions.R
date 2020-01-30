#' Rescale a vector to [0, 1]
#' @description Rescales a vector to [0, 1]
#' @details Rescales a vector to [0, 1].
#' @usage rescale(x)
#' @param \code{x} A vector input
#' @return A vector rescaled to [0, 1]
#' @examples
#' rescale(rnorm(10, 0, 100))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
rescale <- function(x){

  min_x <- min(x)
  range_x <- max(x) - min_x
  scaled <- (x - min_x)/range_x

  return(scaled)

}


#' Scale a matrix to the unit hypercube
#' @description Scales a matrix to \deqn{[0, 1]^d}
#' @details This function scales a matrix to the unit hypercube by subtracting off
#' the minimum value in each column from each value and dividing by the range of
#' that column.
#' @usage scaleX(x)
#' @param \code{x} An \code{n x d} matrix
#' @return A matrix rescaled to \deqn{[0, 1]^d} with attributes
#' \code{"scaled:minimum"} and \code{"scaled:range"}
#' @examples
#' scaleX(matrix(runif(20, 0, 10), nrow = 10, ncol = 2))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
scaleX <- function(x){

  x_scaled <- apply(x, 2, rescale)
  range_x <- apply(x, 2, range)
  x_scaled <- structure(x_scaled,
                        'scaled:minimum' = range_x[1,],
                        'scaled:range' = range_x[2,] - range_x[1,])
  return(x_scaled)
}

#' Rescale a vector
#' @description Rescales a vector to the same scale as a vector that has already
#' been scaled to [0,1]
#' @details Rescales a vector to the same scale as a vector that has already
#' been scaled to [0,1]
#' @usage rescaleXPred(newdata, x)
#' @param \code{newdata} A vector input
#' @param \code{x} A vector input that is scaled to
#' @return A vector rescaled to [0, 1]
#' @examples
#' rescaleXPred(rnorm(10, 0, 100))
#' @author Casey Davis (\email{cbdavis33@@gmail.com})
rescaleXPred <- function(i, newdata, x){
  ## TODO: check the documentation and usage for this
  min_x <- attributes(x)$`scaled:minimum`[i]
  range_x <- attributes(x)$`scaled:range`[i]
  scaled <- (newdata[, i] - min_x)/range_x

  return(scaled)

}
