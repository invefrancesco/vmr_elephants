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
sim.burst <- function(n, mod, x = NULL, burn = 500) {
  # --- Latent variable z ---
  z <- numeric(n)
  if (is.null(x)) {
    z <- sample(1:mod$K, n, replace = TRUE, mod$prob)
  } else {
    linpred <- exp(x %*% t(mod$beta.array)) # (n x K-1)
    den <- 1 + rowSums(linpred) # denominator
    pis <- cbind(1, linpred) / den # probabilities
    z <- apply(pis, 1, \(i) sample(1:mod$K, size = 1, prob = i))
  }

  # --- Init time varying parameters ---
  kappa.t <- matrix(ncol = mod$K, nrow = n + burn)
  kappa.t[1:mod$h, ] <- matrix(mod$kappa, nrow = mod$h, ncol = mod$K, byrow = TRUE)

  mu.t <- matrix(ncol = mod$K, nrow = n + burn)
  mu.t[1:mod$h, ] <- matrix(mod$mu, nrow = mod$h, ncol = mod$K, byrow = TRUE)

  # --- Sampling for the first h time steps ---
  y <- numeric(n + burn)
  for (i in 1:mod$h) {
    y[i] <- circular::rvonmises(1, mu = mod$mu[z[i]], kappa = mod$kappa[z[i]])
  }

  # --- Sampling for the time steps (h+1):burn ---
  for (t in (mod$h + 1):burn) {
    for (k in 1:mod$K) {
      foo <- y[(t - 1):(t - mod$h)]
      foo <- as.numeric(mod$arcoef[, k] %*% sin(foo - mod$mu[k]))
      kappa.t[t, k] <- sqrt(mod$kappa[k]^2 + foo^2)
      mu.t[t, k] <- mod$mu[k] + atan(foo / mod$kappa[k])
    }
    y[t] <- circular::rvonmises(1, mu = mu.t[t, 1], kappa = kappa.t[t, 1])
  }

  # --- Sampling for the time steps (burn + 1):n ---
  for (t in (burn + 1):(burn + n)) {
    for (k in 1:mod$K) {
      foo <- y[(t - 1):(t - mod$h)]
      foo <- as.numeric(mod$arcoef[, k] %*% sin(foo - mod$mu[k]))
      kappa.t[t, k] <- sqrt(mod$kappa[k]^2 + foo^2)
      mu.t[t, k] <- mod$mu[k] + atan(foo / mod$kappa[k])
    }
    y[t] <- circular::rvonmises(1, mu = mu.t[t, z[t - burn]], kappa = kappa.t[t, z[t - burn]])
  }

  y <- y[(burn + 1):length(y)]
  list(y = as.numeric(y), z = z)
}

# --- sim.data ---
#' Generate a dataset from a Circular MAR model
#'
sim.data <- function(dat, mod, formula, ...) {
  dat$xmat <- if (!is.null(formula)) {
    model.matrix(formula, data = dat)
  } else {
    NULL
  }

  dplyr::reframe(dat,
    {
      bind_cols(
        pick(everything()),
        suppressWarnings(
          dplyr::as_tibble(sim.burst(n = n(), mod = mod, x = xmat))
        )
      )
    },
    .by = c(...)
  ) %>%
    mutate(t = row_number(), .by = c(...))
}
