check_seed <- function(seed){

  if(!is.numeric(seed)){
    warning("Seed needs to be numeric. Randomly generating seed.")
    seed <- sample.int(.Machine$integer.max, 1)
  }else if(is.infinite(seed)){
    warning("Seed cannot be Inf. Randomly generating seed.")
    seed <- sample.int(.Machine$integer.max, 1)
  }else{
    seed <- as.integer(seed)
  }
  return(seed)
}


check_valid_cor_mat <- function(x){

  if(!is.matrix(x)) return(FALSE) # check that it's a matrix
  if(!is.numeric(x)) return(FALSE) # check that it's numeric
  if(!(nrow(x) == ncol(x))) return(FALSE) # check square
  if(isFALSE(all(diag(x) == 1))) return(FALSE) # check 1's on diagonal
  if(!(sum(x == t(x)) == (nrow(x)^2))) return(FALSE) # check symmetric

  # check positive semi-definite
  eigenvalues <- eigen(x, only.values = TRUE)$values
  eigenvalues[abs(eigenvalues) < 1e-8] <- 0
  if (any(eigenvalues < 0)) {
    return(FALSE)
  }
  return(TRUE)

}

get_sampler_args_stan <- function(x){

  list(algorithm = "NUTS",
       iter = x@stan_args[[1]]$iter,
       warmup = x@stan_args[[1]]$warmup,
       thin = x@stan_args[[1]]$thin,
       seed = sapply(x@stan_args, function(z) z$seed),
       control = attr(x@sim$samples[[1]], "args")$control)

}

check_valid_x_pred_matrix <- function(x_pred, x){

  if(!is.matrix(x_pred)){
    stop("'x_pred' should be a matrix.")
  }
  if(ncol(x_pred) != ncol(x)){
    stop("'x' and 'x_pred' must have the same number of columns.")
  }
  if(any(is.na(x_pred))){
    stop("'x_pred' should not have any NA values.")
  }

}

coda_summary <- function(samples, quantiles){

  n_eff <- as.matrix(floor(coda::effectiveSize(samples)))
  colnames(n_eff) <- "n_eff"
  R_hats <- as.matrix(coda::gelman.diag(samples)$psrf[, "Point est."])
  colnames(R_hats) <- "Rhat"
  long_sum <- summary(samples, quantiles)

  cbind(long_sum$statistics, long_sum$quantiles,
        n_eff, R_hats)

}

check_valid_truth <- function(object, truth){

  if(!is.numeric(truth)){
    stop("'truth' must be numeric.")
  }
  if(length(truth) != nrow(object@preds$x)){
    stop(strwrap(prefix = " ", initial = "",
                 "The length of 'truth' is not the same as the length
                           of the predictions."))
  }
  if(any(is.na(truth))){
    stop("'truth' must not contain any NA values.")
  }

}

check_prob <- function(x){

  if(!is.numeric(x)){
    stop("'x' must be numeric.")
  }
  if(length(x) != 1){
    stop("'x' must have length 1.")
  }
  if(!is.null(x) && (x <= 0 || x >= 1)){
    stop(strwrap(prefix = " ", initial = "",
                 "x must be in the interval (0, 1)."))
  }

}
