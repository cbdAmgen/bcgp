rm(list = ls())
cat("\014")

# library(bcgp)

# Need to put these in formal tests in tests/ directory

###########################################################
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

plot(data_sim)                       # should work fine
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should work fine with the decomposition
                                     #      plotted
plot(data_sim, raster = TRUE)        # should work fine. Should add a message
                                     #      saying 'raster' is ignored for d = 1
############################################################

rm(list = ls())
cat("\014")

composite <- FALSE
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

plot(data_sim)                       # should work fine
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should return custom error message
plot(data_sim, raster = TRUE)        # should work fine. Should add a message
                                     #      saying 'raster' is ignored for d = 1

############################################################

rm(list = ls())
cat("\014")

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
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE)

plot(data_sim)                       # should work fine
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should work fine with the decomposition
                                     #      plotted
plot(data_sim, raster = TRUE)        # should work fine. Should add a message
                                     #      saying 'raster' is ignored for d = 1

############################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- TRUE
noise <- FALSE
d <- 1
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE)

plot(data_sim)                       # should work fine with a custom warning
plot(data_sim, process = "y")        # should work fine
plot(data_sim, process = "variance") # should give warning about not plotting
                                     #      the variance if stationary
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list (variance
                                     #      plot should be blank)
plot(data_sim, decomposition = TRUE) # should work fine with the decomposition
                                     #      plotted. Custom message for no
                                     #      variance plot
plot(data_sim, process = "y",        # should work fine. Should add a message
     raster = TRUE)                  #      saying 'raster' is ignored for d = 1

############################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- TRUE
d <- 1
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE)

plot(data_sim)                       # should work fine
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should work fine with the decomposition
                                     #      plotted
plot(data_sim, raster = TRUE)        # should work fine. Should add a message
                                     #      saying 'raster' is ignored for d = 1

############################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- TRUE
d <- 1
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

parameters$w <- 1.5
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE) # should throw custom error

parameters$w <- 0.75
parameters$rhoG <- -0.5
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE) # should throw custom error

parameters$rhoG <- 0.5
parameters$rhoL <- 0.7
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE) # should throw custom error

parameters$rhoL <- 0.3
parameters$sigma2V <- -0.4
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE) # should throw custom error

parameters$sigma2V <- 0.1
parameters$rhoV <- 1.5
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE) # should throw custom error

parameters$rhoV <- 0.5
parameters$sigma2eps <- -0.5
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed,
                     randomX1D = TRUE) # should throw custom error

#####################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 2
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 20, n_test = 100, parameters = parameters, seed = seed,
                     grid_test = TRUE)

plot(data_sim)                       # should work fine
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
#      as components of the list
plot(data_sim, decomposition = TRUE) # should give warning and then just plot
                                     #      the overall process and var. process
plot(data_sim, raster = TRUE)        # should work fine


###################################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 2
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 20, n_test = 100, parameters = parameters, seed = seed,
                     grid_test = FALSE)

plot(data_sim)                       # should work fine
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should give warning and then just plot
                                     #      the overall process and var. process
plot(data_sim, raster = TRUE)        # should work fine, but raster is igmored
                                     #      when grid_test = FALSE (maybe have
                                     #      message for that)

###################################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 2
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 20, n_test = 40, parameters = parameters, seed = seed,
                     grid_test = TRUE)

plot(data_sim)                       # should work fine. n_test is ignored.
                                     #     Should give message about that
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should give warning and then just plot
                                     #      the overall process and var. process
plot(data_sim, raster = TRUE)        # should work fine,


###################################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 2
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 20, n_test = 40, parameters = parameters, seed = seed,
                     grid_test = TRUE, grid_test_size = 3)

plot(data_sim)                       # should work fine. n_test is ignored.
                                     #     Should give message about that
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should give warning and then just plot
                                     #      the overall process and var. process
plot(data_sim, raster = TRUE)        # should work fine, but might look funny
                                     #      with the small grid_test_size

###################################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 2
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 20, n_test = 40, parameters = parameters, seed = seed,
                     grid_test = TRUE, grid_test_size = 30)

plot(data_sim)                       # should work fine. n_test is ignored.
                                     #     Should give message about that
plot(data_sim, process = "y")        # should make just the process plot
plot(data_sim, process = "variance") # should make just the variance plot
z <- plot(data_sim, print = FALSE)   # should return a list with the two plots
                                     #      as components of the list
plot(data_sim, decomposition = TRUE) # should give warning and then just plot
                                     #      the overall process and var. process
plot(data_sim, raster = TRUE)        # should work fine

###################################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 6
seed <- 12345

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 20, n_test = 40, parameters = parameters, seed = seed,
                     grid_test = TRUE, grid_test_size = 3)

plot(data_sim)  # Should give custom error message about plotting not supported


data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 20, n_test = 40, parameters = parameters, seed = seed,
                     grid_test = TRUE) # Give warning about grid_test_size

###########################################################

rm(list = ls())
cat("\014")

composite <- TRUE
stationary <- FALSE
noise <- FALSE
d <- 1
seed <- "xyz"

parameters <- create_parameter_list(composite = composite,
                                    stationary = stationary,
                                    noise = noise, d = d)

# Should give warning about seed
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed)

plot(data_sim)                       # should work fine

seed <- -5
data_sim <- bcgpsims(composite = composite, stationary = stationary,
                     noise = noise, d = d,
                     n = 12, n_test = 50, parameters = parameters, seed = seed)

plot(data_sim)                       # should work fine
