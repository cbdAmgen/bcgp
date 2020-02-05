.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}

.onAttach <- function(...) {
  bcgpLib <- dirname(system.file(package = "bcgp"))
  # pkgdesc <- packageDescription("bcgp", lib.loc = bcgpLib)
  pkgdesc <- utils::packageDescription("bcgp")
  packageStartupMessage(paste("bcgp (Version ", pkgdesc$Version,")", sep = ""))
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM call\n",
                        "options(mc.cores = parallel::detectCores()).")
}
