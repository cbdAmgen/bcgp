central_intervals <- function(object, prob){

  stopifnot(is.matrix(object))
  if (length(prob) != 1L || prob <= 0 || prob >= 1)
    stop("'prob' should be a single number greater than 0 and less than 1.",
         call. = FALSE)
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  labs <- paste0(100 * probs, "%")
  out <- t(apply(object, 2, quantile, probs = probs))
  dimnames(out) <- list(colnames(object), c("lower", "upper"))
  attr(out, "Probability") <- prob

  return(out)

}

hpd_intervals <- function(object, prob){

  stopifnot(is.matrix(object))
  if (length(prob) != 1L || prob <= 0 || prob >= 1)
    stop("'prob' should be a single number greater than 0 and less than 1.",
         call. = FALSE)

  return(coda::HPDinterval(coda::as.mcmc(object), prob = prob))

}
