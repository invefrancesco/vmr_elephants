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
prob <- c(0.25, 0.75) # Mixing probabilities
kappa <- c(5, 8) # Unconditional concentrations
mu <- c(0, pi) # Unconditional means

# A matrix h * K where h is the autoregressive order and K is the number of
# mixture components
arcoef <- matrix( # Autoregressive parameters
  c(0.1, 0.1),
  ncol = 2, byrow = TRUE
)

set.seed(1234)

burst <- rep(1:100, each = 500)
x <- as.vector(replicate(100, simulate_burst(500, prob, kappa, mu, arcoef)))

fit <- fit_cmar(x, burst, K = 2, h = 1, tol = 1e-4)
summary(fit)
plot(fit$history)

# TODO
# [ ]
