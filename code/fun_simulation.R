#' # data simulation
#'
#' this script implement the functions to simulate data from a circular mar
#' model using von mises distributions.

pacman::p_load(
  circular
)
#' --------------------------------------------------------------------------- #
#'   MAIN SIMULATION FUNCTION
#' --------------------------------------------------------------------------- #
#' --- sim.burst ---
#' Generate a single trajectory from a Circular MAR model
#'
#' @param n number of observations (after burn-in)
#' @param mod list containing model parameters:
#'   - K: number of components
#'   - h: autoregressive order
#'   - mu: baseline mean directions (length K)
#'   - kappa: baseline concentrations (length K)
#'   - arcoef: AR coefficients (matrix h x K)
#'   - prob: mixing probabilities (length K)
#' @param burn length of burn-in period to reach stationarity
#'
#' @return a vector of circular observations
sim.burst <- function(n, mod, burn = 500) {
  n <- n + burn
  # --- Init time varying parameters ---
  kappa.t <- matrix(ncol = mod$K, nrow = n)
  kappa.t[1:mod$h, ] <- matrix(mod$kappa, nrow = mod$h, ncol = mod$K, byrow = TRUE)

  mu.t <- matrix(ncol = mod$K, nrow = n)
  mu.t[1:mod$h, ] <- matrix(mod$mu, nrow = mod$h, ncol = mod$K, byrow = TRUE)

  # --- Latent variable z ---
  z <- sample(1:mod$K, n, replace = TRUE, mod$prob)

  # --- Initial sampling for the first h time steps ---
  x <- numeric(n)
  for (i in 1:mod$h) {
    x[i] <- circular::rvonmises(1, mu = mod$mu[z[i]], kappa = mod$kappa[z[i]])
  }

  # --- Sampling for the time steps h:n ---
  for (t in (mod$h + 1):n) {
    for (k in 1:mod$K) {
      foo <- x[(t - 1):(t - mod$h)]
      foo <- as.numeric(mod$arcoef[, k] %*% sin(foo - mod$mu[k]))
      kappa.t[t, k] <- sqrt(mod$kappa[k]^2 + foo^2)
      mu.t[t, k] <- mod$mu[k] + atan(foo / mod$kappa[k])
    }
    x[t] <- circular::rvonmises(1, mu = mu.t[t, z[t]], kappa = kappa.t[t, z[t]])
  }
  x <- x[(burn + 1):length(x)]
  z <- z[(burn + 1):length(z)]
  list(x = x, z = z)
}

n <- rpois(5, 10)

# --- sim.data ---
#' Generate a dataset from a Circular MAR model
#'
sim.data <- function(n, mod) {
  dat <- lapply(seq_along(n), function(b) {
    sim <- sim.burst(n = n[b], mod = mod)
    data.frame(
      burst = rep(b, n[b]),
      x = sim$x,
      z = sim$z
    )
  })
  do.call(rbind, dat)
}
