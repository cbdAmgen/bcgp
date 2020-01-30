validate_parameter_list <- function(parameters, composite, stationary, d){

  if(isTRUE(composite)){
    if(isFALSE(stationary)){
      validate_params_comp_ns(parameters, d)
    }else{ # composite == TRUE, stationary == TRUE
      validate_params_comp_s(parameters, d)
    }
  }else{
    if(isFALSE(stationary)){
      validate_params_noncomp_ns(parameters, d)
    }else{ # composite == FALSE, stationary == TRUE
      validate_params_noncomp_s(parameters, d)
    }
  }

  return(NULL)
}

validate_params_comp_ns <- function(parameters, d){

  with(parameters, stopifnot(length(rhoG) == d,
                             length(rhoL) == d,
                             length(rhoV) == d,
                             is.numeric(beta0),
                             (w >= 0 & w <= 1),
                             all(0 <= rhoG & rhoG <= 1),
                             all(0 <= rhoL & rhoL <= rhoG),
                             is.numeric(muV),
                             (sigma2V >= 0),
                             all(0 <= rhoV & rhoV <= 1),
                             sigma2eps >= 0))
  return(NULL)
}

validate_params_comp_s <- function(parameters, d){

  with(parameters, stopifnot(length(rhoG) == d,
                             length(rhoL) == d,
                             is.numeric(beta0),
                             (w >= 0 & w <= 1),
                             all(0 <= rhoG & rhoG <= 1),
                             all(0 <= rhoL & rhoL <= rhoG),
                             (sigma2 >= 0),
                             sigma2eps >= 0))
  return(NULL)
}

validate_params_noncomp_ns <- function(parameters, d){

  with(parameters, stopifnot(length(rho) == d,
                             length(rhoV) == d,
                             is.numeric(beta0),
                             all(0 <= rho & rho <= 1),
                             is.numeric(muV),
                             (sigma2V >= 0),
                             all(0 <= rhoV & rhoV <= 1),
                             sigma2eps >= 0))
  return(NULL)
}

validate_params_noncomp_s <- function(parameters, d){

  with(parameters, stopifnot(length(rho) == d,
                             is.numeric(beta0),
                             all(0 <= rho & rho <= 1),
                             (sigma2 >= 0),
                             sigma2eps >= 0))
  return(NULL)
}
