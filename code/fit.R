pacman::p_load(
  circular,
  tidyverse,
  sf
)
source("code/fun_estimation.R")

fit <- function(nSim, K, h) {
  load(paste0("data/sim", nSim, ".RData"))
  fit <- cmar.ms(dat$x, dat$burst, K, h)
  save(fit, file = paste0("data/fit", nSim, ".RData"))
}

fit(6, 4, 1)
