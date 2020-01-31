validate_bcgpmodel_inputs <- function(x, y, composite, stationary, noise){

  validate_data(x, y)
  validate_logical(composite)
  validate_logical(stationary)
  validate_logical(noise)

}

validate_data <- function(x, y){

  validate_x(x)
  validate_y(y)
  validate_xy(x, y)

}

validate_x <- function(x){

  if(!is.matrix(x)) stop("'x' must be a matrix.")
  if(!is.numeric(x)) stop("'x' must be numeric.")
  if(any(is.na(x))) stop("'x' must not have any NA values.")

}

validate_y <- function(y){

  if(!is.numeric(y)) stop("'y' must be numeric.")
  if(any(is.na(y))) stop("'y' must not have any NA values.")

}

validate_xy <- function(x, y){

  if(nrow(x) != length(y))
    stop("'x' must have the same number of rows as length of 'y'.")

}

validate_logical <- function(x){

  if(!(isTRUE(x) || isFALSE(x)))
    stop(strwrap(prefix = " ", initial = "",
                 paste0("'", deparse(substitute(x)), "'", " must be either
                        'TRUE' or 'FALSE'.")))

}
