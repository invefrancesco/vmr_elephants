#' # circular mixture autoregressive (cmar) model
#'
#' this is the main script.

pacman::p_load(
  circular
)
source("code/fun_estimation.R")

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
  c(
    0.5, 0.5,
    0.2, 0.2
  ),
  ncol = 2, byrow = TRUE
)

set.seed(1234)

burst <- rep(1:5, each = 100)
x <- as.vector(replicate(5, simulate_burst(100, prob, kappa, mu, arcoef)))

fit <- fit_cmar(x, burst, K = 2, h = 1, tol = 1e-6)
summary(fit)
plot(fit$history)

# TODO
# [ ]
