rm(list = ls())
cat("\014")

library(bcgp)

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 1
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed)

plot(data_sim)

my_model_comp_ns <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
                              composite = TRUE, stationary = FALSE,
                              noise = FALSE)

## If you want to change the priors
# my_model@priors$w$alpha <- 2
# my_model@priors$w

my_fit_comp_ns <- bcgp_sampling(my_model_comp_ns, algorithm = "NUTS",
                                scaled = TRUE, chains = 4L,
                                cores = 4L, init = "random", burnin = 1000,
                                n_mcmc = 2000,
                                control = list(adapt_delta = 0.99,
                                               max_treedepth = 13))

summary(my_fit_comp_ns)
print(my_fit_comp_ns)
my_fit_comp_ns

##################################################

my_model_comp_s <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
                             composite = TRUE, stationary = TRUE, noise = FALSE)

my_fit_comp_s <- bcgp_sampling(my_model_comp_s, algorithm = "NUTS",
                               scaled = TRUE, chains = 4L,
                               cores = 4L, init = "random", burnin = 1000,
                               n_mcmc = 2000,
                               control = list(adapt_delta = 0.99,
                                              max_treedepth = 13))

summary(my_fit_comp_s)
print(my_fit_comp_s)
my_fit_comp_s

####################################################

my_model_noncomp_s <- bcgpmodel(x = data_sim@training$x,
                                y = data_sim@training$y,
                                composite = FALSE, stationary = TRUE,
                                noise = FALSE)

my_fit_noncomp_s <- bcgp_sampling(my_model_noncomp_s, algorithm = "NUTS",
                                  scaled = TRUE, chains = 4L,
                                  cores = 4L, init = "random", burnin = 1000,
                                  n_mcmc = 2000,
                                  control = list(adapt_delta = 0.99,
                                                 max_treedepth = 13))

summary(my_fit_noncomp_s)
print(my_fit_noncomp_s)
my_fit_noncomp_s
####################################################

my_model_noncomp_ns <- bcgpmodel(x = data_sim@training$x,
                                 y = data_sim@training$y,
                                 composite = FALSE, stationary = FALSE,
                                 noise = FALSE)

my_fit_noncomp_ns <- bcgp_sampling(my_model_noncomp_ns, algorithm = "NUTS",
                                   scaled = TRUE, chains = 4L,
                                   cores = 4L, init = "random", burnin = 1000,
                                   n_mcmc = 2000,
                                   control = list(adapt_delta = 0.99,
                                                  max_treedepth = 13))

summary(my_fit_noncomp_ns)
print(my_fit_noncomp_ns)
my_fit_noncomp_ns

print(my_fit_noncomp_ns, pars = "V")
