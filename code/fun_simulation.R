#' # data simulation
#'
#' this script implement the functions to simulate data from a circular mar
#' model using von mises distributions.

pacman::p_load(
  circular
)

# --------------------------------------------------------------------------- #
#   AUXILIARY FUNCTIONS
# --------------------------------------------------------------------------- #
# --- get_tau ---
#'
#' @description
#' Gives a matrix where each cell contains the circular deviation at lag h for
#' component k.
#'
#' @param t the current time
#' @param K the number of mixture components
#' @param h the autoregressive order
#' @param mu the non-conditional mean

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

# --- get_upsilon ---
#'
#' @description
#' Gives a matrix where each cell contains the circular deviation at lag h for
#' component k scaled by concentration kappa.
#'
#' @param t the current time
#' @param K the number of mixture components
#' @param h the autoregressive order
#' @param mu the non-conditional mean
#'
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
# --- simulate_burst ---
#'
#'
simulate_burst <- function(n, prob, kappa, mu, arcoef) {
  K <- length(prob) # Number of mixture components
  h <- nrow(arcoef) # Autoregressive order

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

  # --- Sampling for the time steps 2:n ---
  for (t in 2:n) {
    tau <- get_tau(t, K, h, x, mu)
    kappa.t[t, ] <- sqrt(kappa^2 + colSums(arcoef * tau)^2)

    upsilon <- get_upsilon(t, K, h, x, mu, kappa)
    mu.t[t, ] <- mu + atan(colSums(arcoef * upsilon))

    x[t] <- circular::rvonmises(1, mu = mu.t[t, z[t]], kappa = kappa.t[t, z[t]])
  }

  x
}
