#' # circular mixture autoregressive (cmar) model
#'
#' this is the main script.

pacman::p_load(
  circular
)
source("code/fun_estimation.R")
source("code/fun_simulation.R")

# --------------------------------------------------------------------------- #
#   TRIAL PARAMETERS AND SIMULATION
# --------------------------------------------------------------------------- #

# --- Trial parameters ---
K <- 2
h <- 3
mod <- list(
  K = K,
  h = h,
  mu = c(0, pi),
  kappa = c(10, 12),
  prob = c(0.4, 0.6),
  arcoef = matrix(0.1, nrow = h, ncol = K)
)

# 3. Ora puoi lanciare la simulazione
set.seed(1234)
x <- sim.burst(n = 10, mod = mod)

burst <- rep(1:5, each = 10)
x <- as.vector(replicate(5, sim.burst(n = 10, mod = mod, burn = 500)))


fit <- fit_cmar(x, burst, K = 2, h = 1, tol = 1e-4)
summary(fit)
plot(fit$history)

# TODO
# [ ] Salvare la componente reale in sim.burst$z
