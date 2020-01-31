rm(list = ls())
cat("\014")

# library(bcgp)

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

my_model <- bcgpmodel(x = data_sim@training$x, y = data_sim@training$y,
                      composite = TRUE, stationary = FALSE, noise = FALSE)

## If you want to change the priors
# my_model@priors$w$alpha <- 2
# my_model@priors$w

