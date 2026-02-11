#' # Simulation
#'
#' This file contains the code to simulate the data

rm(list = ls())
pacman::p_load(
  circular
)

# --------------------------------------------------------------------------- #
#   AUXILIARY FUNCTIONS
# --------------------------------------------------------------------------- #
# COMMENT: explain the matrix tau
get_tau <- function(t, K, h, x, mu) {
  tau <- matrix(ncol = K, nrow = h)
  for (g in 1:h) {
    ifelse(
      t > g,
      tau[g, ] <- sin(x[t - g] - mu),
      tau[g, ] <- 0
    )
  }
  tau
}

# COMMENT: explain the matrix upsilon
get_upsilon <- function(t, K, h, x, mu, kappa) {
  upsilon <- matrix(ncol = K, nrow = h)
  for (g in 1:h) {
    ifelse(
      t > g,
      upsilon[g, ] <- sin(x[t - g] - mu) / kappa,
      upsilon[g, ] <- 0
    )
  }
  upsilon
}

# --------------------------------------------------------------------------- #
#   MAIN SIMULATION FUNCTION
# --------------------------------------------------------------------------- #
simulate_burst <- function(n, prob, kappa, mu, phi) {
  K <- length(prob) # Number of mixture components
  h <- nrow(phi) # Autoregressive order

  # --- Init time varying parameters ---
  kappa.t <- matrix(ncol = K, nrow = n)
  kappa.t[1, ] <- kappa

  mu.t <- matrix(ncol = K, nrow = n)
  mu.t[1, ] <- mu

  # --- Latent variable z ---
  z <- sample(1:K, n, replace = TRUE, prob)

  # --- Initial sampling for the first time step (t=1) ---
  x <- numeric(n)
  x[1] <- circular::rvonmises(1, mu = mu[z[1]], kappa = kappa[z[1]])

  # --- Title here ---
  for (t in 2:n) {
    tau <- get_tau((t - 1), K, h, x, mu)
    kappa.t[t, ] <- sqrt(kappa^2 + colSums(phi * tau)^2)

    upsilon <- get_upsilon((t - 1), K, h, x, mu, kappa)
    mu.t[t, ] <- mu + atan(colSums(phi * upsilon))

    x[t] <- circular::rvonmises(1, mu = mu.t[t, z[t]], kappa = kappa.t[t, z[t]])
  }

  x
}

# --------------------------------------------------------------------------- #
#   TRIAL PARAMETERS AND SIMULATION
# --------------------------------------------------------------------------- #

n <- 100 # Total number of observations to simulate
prob <- c(0.25, 0.75) # Mixing probabilities for the latent states (Z)

kappa <- c(8, 10) # Unconditional concentrations
mu <- c(0, pi) # Unconditional means

# A matrix h * K where h is the autoregressive order and K is the number of
# mixture components
phi <- matrix( # Autoregressive parameters
  c(
    0.05, 0.05,
    0.025, 0.025,
    0.01, 0.01
  ),
  ncol = 2, byrow = TRUE
)

set.seed(1234)
simulate_burst(100, prob, kappa, mu, phi) # Trial simulation
