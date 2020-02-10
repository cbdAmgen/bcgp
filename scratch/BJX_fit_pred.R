rm(list = ls())
cat("\014")

library(bcgp)

data(BJX)

x <- BJX$x
y <- BJX$y

x_pred <- BJX$x_test
y_test <- BJX$y_test

model_c_ns_with_pred <- bcgpmodel(x = x, y = y,
                                  composite = TRUE,
                                  stationary = FALSE,
                                  noise = FALSE)

model_c_ns_with_pred@priors$sigma2V$alpha <- 3
model_c_ns_with_pred@priors$sigma2V$beta <- 3

fit_c_ns_with_pred <-
  bcgp_sampling_and_prediction(model_c_ns_with_pred,
                               x_pred = x_pred,
                               algorithm = "NUTS",
                               scaled = FALSE,
                               chains = 4L,
                               cores = 4L,
                               init = "random",
                               burnin = 1000,
                               n_mcmc = 2000,
                               control = list(adapt_delta = 0.9925,
                                              max_treedepth = 15))

print(fit_c_ns_with_pred, digits_summary = 3,
      pars = c("beta0", "w", "rhoG", "rhoL", "sigma2eps", "muV", "rhoV",
               "sigma2V"))

print(fit_c_ns_with_pred, digits_summary = 3)

summary(fit_c_ns_with_pred)
