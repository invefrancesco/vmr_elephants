#' # EM Estimation
#'
#' This file contains the functions needed to compute the estimations

rm(list = ls())
pacman::p_load(
  circular
)

# --------------------------------------------------------------------------- #
#   AUXILIARY FUNCTIONS
# --------------------------------------------------------------------------- #
# Tau Matrix: Captures the raw circular deviations between the current
# observation and its past lags.
get_tau <- function(t, K, h, x, mu) {
  # K: number of mixture component in the model
  # h: maximum lag order (memory depth) of the autoregressive process
  tau <- matrix(ncol = K, nrow = h)
  for (g in 1:h) {
    ifelse(
      t > g, # Checks if enough historical data exists for the specific lag 'g'
      tau[g, ] <- sin(x[t - g] - mu),
      tau[g, ] <- 0
    )
  }
  tau
}

# Upsilon Matrix: Computes the past circular deviations scaled by the
# concentration parameter (kappa). It modulates the impact of historical lags.
get_upsilon <- function(t, K, h, x, mu, kappa) {
  upsilon <- matrix(ncol = K, nrow = h)
  for (g in 1:h) {
    ifelse(
      t > g, # Checks if enough historical data exists for the specific lag 'g'
      upsilon[g, ] <- sin(x[t - g] - mu) / kappa,
      upsilon[g, ] <- 0
    )
  }
  upsilon
}

# --------------------------------------------------------------------------- #
#   LIKELIHOOD FUNCTION
# --------------------------------------------------------------------------- #


loglik <- function(x, burst, prob, mu, kappa, arcoef, K, h) {
  # Initialize list to compute dynamic parameters
  kappa.t <- list()
  mu.t <- list()

  # Compute dynamic parameters independently for each trajectory
  for (b in unique(burst)) {
    x.b <- x[burst == b]

    # Pre-allocate temporary matrices for the current burst 'b'
    kappa.tmp <- matrix(nrow = length(x.b), ncol = K)
    mu.tmp <- matrix(nrow = length(x.b), ncol = K)

    for (t in seq_along(x.b)) {
      tau <- get_tau(t, K, h, x.b, mu)
      upsilon <- get_upsilon((t - 1), K, h, x.b, mu, kappa)

      # Conditional update of concentration and mean for time t
      kappa.tmp[t, ] <- sqrt(kappa^2 + colSums(arcoef * tau)^2)
      mu.tmp[t, ] <- mu + atan(colSums(arcoef * upsilon))
    }

    # Store the burst-specific matrices into the main lists
    kappa.t[[b]] <- kappa.tmp
    mu.t[[b]] <- mu.tmp
  }

  # Flatten the lists: row-binds the separated bursts into continuous matrices
  # to align perfectly with the original data vector 'x'.
  kappa.t <- do.call(rbind, kappa.t)
  mu.t <- do.call(rbind, mu.t)

  # Compute the likelihood matrix (l).
  # Each cell (i, k) is the probability of observing x[i] in state k,
  # weighted by the state's prior probability (prob[k]).
  l <- matrix(ncol = K, nrow = length(x))
  for (k in seq_len(K)) {
    for (i in seq_along(x)) {
      l[i, k] <- prob[k] * circular::dvonmises(x[i], mu.t[i, k], kappa.t[i, k])
    }
  }

  # Return the likelihood matrix and the total log-likelihood
  list(
    l,
    sum(log(rowSums(l)))
  )
}

# --------------------------------------------------------------------------- #
#   MAIN FUNCTION (EM ALGORITHM)
# --------------------------------------------------------------------------- #

prob <- rexp(2)
prob <- init_prob / sum(init_prob)
mu <- runif(2, 0, 2 * pi)
kappa <- rexp(2)
arcoef <- runif(2 * 2)
x <- dat[, 2]
burst <- dat[, 1]

K <- length(prob)
h <- length(arcoef)

l <- vector()
l <- -Inf
l[2] <- loglik(x, burst, prob, mu, kappa, arcoef, K, h)[[2]]

# E step
f <- loglik(x, burst, prob, mu, kappa, arcoef, K, h)[[1]]
w <- f / rowSums(f)

# M step
