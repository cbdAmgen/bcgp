bcgp_stan <- function(x, scaled, chains, cores, iter, warmup, thin, control,
                      ...){

  # dots <- list(...)
  # if("thin" %in% names(dots) && dots$thin != 1){
  #   message(strwrap(prefix = " ", initial = "",
  #                   "There's no reason to thin these models."))
  #
  # }

  if(isTRUE(x@composite)){
    if(isFALSE(x@stationary)){
      ## composite, non- stationary
      out <- bcgp_stan_comp_ns(x, scaled, chains, cores, iter, warmup, thin,
                               control, ...)
    }else{
      ## composite, stationary
      out <- bcgp_stan_comp_s(x, scaled, chains, cores, iter, warmup, thin,
                              control, ...)
    }
  }else{
    if(isFALSE(x@stationary)){
      ## non-composite, non- stationary
      out <- bcgp_stan_noncomp_ns(x, scaled, chains, cores, iter, warmup, thin,
                                  control, ...)
    }else{
      ## non-composite, stationary
      out <- bcgp_stan_noncomp_s(x, scaled, chains, cores, iter, warmup, thin,
                                 control, ...)
    }
  }
  return(out)
}

bcgp_stan_noncomp_s <- function(x, scaled, chains, cores, iter, warmup, thin,
                                control, ...){


  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    x_in <- data$scaled$x
    y_in <- data$scaled$y
  }else{
    x_in <- x@data$x
    y_in <- x@data$y
  }

  d <- ncol(x_in)
  n <- length(y_in)

  stan_data <- list(x = x_in,
                    y = as.vector(y_in),
                    n = n,
                    d = d,
                    rho_alpha = array(x@priors$rho$alpha, dim = d),
                    rho_beta = array(x@priors$rho$beta, dim = d),
                    sigma2_alpha = x@priors$sigma2$alpha,
                    sigma2_beta = x@priors$sigma2$beta,
                    sigma2eps_alpha = x@priors$sigma2eps$alpha,
                    sigma2eps_beta = x@priors$sigma2eps$beta)


  # dots <- list(...)
  # args <- prepare_sampling_args(x, dots)



  stan_fit <- rstan::sampling(stanmodels$stanNonCompS, data = stan_data,
                              pars = c("beta0", "rho", "sigma2eps", "sigma2"),
                              chains = chains, cores = cores, iter = iter,
                              warmup = warmup, thin = thin, control = control,
                              ...)
  # stan_fit <- do.call(rstan::sampling, list(object = stanmodels$stanNonCompS, data = stan_data,
  #                    pars = c("beta0", "rho", "sigma2eps", "sigma2"), thin = dots$thin,
  #                    chains))
  sampler_args <- get_sampler_args_stan(stan_fit)

  # model_name = out1$model_name,
  # data = out1$data,
  # stationary = object@stationary,
  # composite = object@composite,
  # noise = object@noise,
  # scaled = scaled,
  # chains = chains,
  # priors = object@priors,
  # distributions = object@distributions,
  # init = out1$inits,
  # model_pars = out1$model_pars,
  # par_dims = out1$par_dims,
  # sim = out1$sim,
  # algorithm = algorithm,
  # sampler_args = out1$sampler_args,
  # date = date())

  out <- list(data = data,
              model_name = "noncomposite_stationary",
              init = stan_fit@inits,
              model_pars = c("beta0", "rho", "sigma2eps", "sigma2"),
              par_dims = list(beta0 = numeric(0),
                              rho = d,
                              sigma2eps = numeric(0),
                              sigma2 = numeric(0)),
              sims = rstan::As.mcmc.list(stan_fit),
              sampler_args = sampler_args)

  return(out)
}

bcgp_stan_comp_s <- function(x, scaled, chains, cores, iter, warmup, thin,
                             control, ...){

  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    x_in <- data$scaled$x
    y_in <- data$scaled$y
  }else{
    x_in <- x@data$x
    y_in <- x@data$y
  }

  d <- ncol(x_in)
  n <- length(y_in)

  stan_data <- list(x = x_in,
                    y = as.vector(y_in),
                    n = n,
                    d = d,
                    w_lower = x@priors$w$lower, w_upper = x@priors$w$upper,
                    w_alpha = x@priors$w$alpha, w_beta = x@priors$w$beta,
                    rhoG_alpha = array(x@priors$rhoG$alpha, dim = d),
                    rhoG_beta = array(x@priors$rhoG$beta, dim = d),
                    rhoL_alpha = array(x@priors$rhoL$alpha, dim = d),
                    rhoL_beta = array(x@priors$rhoL$beta, dim = d),
                    sigma2_alpha = x@priors$sigma2$alpha,
                    sigma2_beta = x@priors$sigma2$beta,
                    sigma2eps_alpha = x@priors$sigma2eps$alpha,
                    sigma2eps_beta = x@priors$sigma2eps$beta)

  stan_fit <- rstan::sampling(stanmodels$stanCompS, data = stan_data,
                              pars = c("beta0", "w", "rhoG", "rhoL",
                                       "sigma2eps", "sigma2"),
                              chains = chains, cores = cores, iter = iter,
                              warmup = warmup, thin = thin, control = control,
                              ...)

  sampler_args <- get_sampler_args_stan(stan_fit)

  out <- list(data = data,
              model_name = "composite_stationary",
              init = stan_fit@inits,
              model_pars = c("beta0", "w", "rhoG", "rhoL", "sigma2eps",
                             "sigma2"),
              par_dims = list(beta0 = numeric(0),
                              w = numeric(0),
                              rhoG = d,
                              rhoL = d,
                              sigma2eps = numeric(0),
                              sigma2 = numeric(0)),
              sims = rstan::As.mcmc.list(stan_fit),
              sampler_args = sampler_args)

  return(out)
}

bcgp_stan_comp_ns <- function(x, scaled, chains, cores, iter, warmup, thin,
                              control, ...){

  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    x_in <- data$scaled$x
    y_in <- data$scaled$y
  }else{
    x_in <- x@data$x
    y_in <- x@data$y
  }

  d <- ncol(x_in)
  n <- length(y_in)

  stan_data <- list(x = x_in,
                    y = as.vector(y_in),
                    n = n,
                    d = d,
                    w_lower = x@priors$w$lower,
                    w_upper = x@priors$w$upper,
                    w_alpha = x@priors$w$alpha,
                    w_beta = x@priors$w$beta,
                    rhoG_alpha = array(x@priors$rhoG$alpha, dim = d),
                    rhoG_beta = array(x@priors$rhoG$beta, dim = d),
                    rhoL_alpha = array(x@priors$rhoL$alpha, dim = d),
                    rhoL_beta = array(x@priors$rhoL$beta, dim = d),
                    muV_betaV = x@priors$muV$betaV,
                    muV_sigma2 = x@priors$muV$sigma2,
                    rhoV_alpha = array(x@priors$rhoV$alpha, dim = d),
                    rhoV_beta = array(x@priors$rhoV$beta, dim = d),
                    sigma2V_alpha = x@priors$sigma2V$alpha,
                    sigma2V_beta = x@priors$sigma2V$beta,
                    sigma2eps_alpha = x@priors$sigma2eps$alpha,
                    sigma2eps_beta = x@priors$sigma2eps$beta)

  stan_fit <- rstan::sampling(stanmodels$stanCompNS, data = stan_data,
                              pars = c("beta0", "w", "rhoG", "rhoL", "sigma2eps",
                                       "muV", "sigma2V", "rhoV", "V"),
                              chains = chains, cores = cores, iter = iter,
                              warmup = warmup, thin = thin, control = control,
                              ...)

  sampler_args <- get_sampler_args_stan(stan_fit)

  out <- list(data = data,
              model_name = "composite_nonstationary",
              init = stan_fit@inits,
              model_pars = c("beta0", "w", "rhoG", "rhoL", "sigma2eps",
                             "muV", "sigma2V", "rhoV", "V"),
              par_dims = list(beta0 = numeric(0),
                              w = numeric(0),
                              rhoG = d,
                              rhoL = d,
                              sigma2eps = numeric(0),
                              muV = numeric(0),
                              sigma2V = numeric(0),
                              rhoV = d,
                              V = n),
              sims = rstan::As.mcmc.list(stan_fit),
              sampler_args = sampler_args)

  return(out)
}


bcgp_stan_noncomp_ns <- function(x, scaled, chains, cores, iter, warmup, thin,
                                 control, ...){

  data <- list(raw = list(x = x@data$x, y = x@data$y),
               scaled = list(x = scaleX(x@data$x),
                             y = scale(x@data$y, center = TRUE,
                                       scale = TRUE)))


  if(isTRUE(scaled)){
    x_in <- data$scaled$x
    y_in <- data$scaled$y
  }else{
    x_in <- x@data$x
    y_in <- x@data$y
  }

  d <- ncol(x_in)
  n <- length(y_in)

  stan_data <- list(x = x_in,
                    y = as.vector(y_in),
                    n = n,
                    d = d,
                    rho_alpha = array(x@priors$rho$alpha, dim = d),
                    rho_beta = array(x@priors$rho$beta, dim = d),
                    muV_betaV = x@priors$muV$betaV,
                    muV_sigma2 = x@priors$muV$sigma2,
                    rhoV_alpha = array(x@priors$rhoV$alpha, dim = d),
                    rhoV_beta = array(x@priors$rhoV$beta, dim = d),
                    sigma2V_alpha = x@priors$sigma2V$alpha,
                    sigma2V_beta = x@priors$sigma2V$beta,
                    sigma2eps_alpha = x@priors$sigma2eps$alpha,
                    sigma2eps_beta = x@priors$sigma2eps$beta)


  stan_fit <- rstan::sampling(stanmodels$stanNonCompNS, data = stan_data,
                              pars = c("beta0", "rho", "sigma2eps",
                                       "muV", "sigma2V", "rhoV", "V"),
                              chains = chains, cores = cores, iter = iter,
                              warmup = warmup, thin = thin, control = control,
                              ...)

  sampler_args <- get_sampler_args_stan(stan_fit)

  out <- list(data = data,
              model_name = "noncomposite_nonstationary",
              init = stan_fit@inits,
              model_pars = c("beta0", "rho", "sigma2eps",
                             "muV", "sigma2V", "rhoV", "V"),
              par_dims = list(beta0 = numeric(0),
                              rho = d,
                              sigma2eps = numeric(0),
                              muV = numeric(0),
                              sigma2V = numeric(0),
                              rhoV = d,
                              V = n),
              sims = rstan::As.mcmc.list(stan_fit),
              sampler_args = sampler_args)

  return(out)

}
