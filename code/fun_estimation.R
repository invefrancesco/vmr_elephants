#' # circular mixture autoregressive (cmar) model
#'
#' this script implements the expectation-maximization (em) algorithm
#' to estimate parameters of a circular mar model using von mises distributions.
#' it handles independent bursts of circular time-series data.

pacman::p_load(
  circular
)
source("code/fun_simulation.R")

# --------------------------------------------------------------------------- #
#   LIKELIHOOD FUNCTION
# --------------------------------------------------------------------------- #
loglik <- function(param, prob, x, burst, K, h) {
  # --- Unpack parameters ---
  mu <- param[1:K]
  kappa <- exp(param[(K + 1):(2 * K)])
  arcoef <- matrix(
    param[(2 * K + 1):((2 * K) + h * K)],
    nrow = h, ncol = K, byrow = TRUE
  )

  kappa.t <- list()
  mu.t <- list()

  # --- Dynamic parameters per burst ---
  for (b in unique(burst)) {
    x.b <- x[burst == b]

    kappa.tmp <- matrix(nrow = length(x.b), ncol = K)
    mu.tmp <- matrix(nrow = length(x.b), ncol = K)

    for (t in seq_along(x.b)) {
      tau <- get_tau(t, K, h, x.b, mu)
      upsilon <- get_upsilon(t, K, h, x.b, mu, kappa)

      # --- Conditional updates ---
      kappa.tmp[t, ] <- sqrt(kappa^2 + colSums(arcoef * tau)^2)
      mu.tmp[t, ] <- mu + atan(colSums(arcoef * upsilon))
    }

    kappa.t[[b]] <- kappa.tmp
    mu.t[[b]] <- mu.tmp
  }

  # --- Flatten lists to match original vector length ---
  kappa.t <- do.call(rbind, kappa.t)
  mu.t <- do.call(rbind, mu.t)

  # --- Likelihood matrix (observing x[i] given param[k]) ---
  l <- matrix(ncol = K, nrow = length(x))
  for (k in seq_len(K)) {
    for (i in seq_along(x)) {
      l[i, k] <- prob[k] * circular::dvonmises(x[i], mu.t[i, k], kappa.t[i, k])
    }
  }

  list(
    l,
    sum(log(rowSums(l)))
  )
}

# --------------------------------------------------------------------------- #
#   Q FUNCTION
# --------------------------------------------------------------------------- #
Q <- function(param, x, burst, w, K, h) {
  # --- Unpack parameters ---
  mu <- param[1:K]
  kappa <- exp(param[(K + 1):(2 * K)])
  arcoef <- matrix(
    param[(2 * K + 1):((2 * K) + h * K)],
    nrow = h, ncol = K, byrow = TRUE
  )

  kappa.t <- list()
  mu.t <- list()

  # --- Dynamic parameters per burst ---
  for (b in unique(burst)) {
    x.b <- x[burst == b]

    kappa.tmp <- matrix(nrow = length(x.b), ncol = K)
    mu.tmp <- matrix(nrow = length(x.b), ncol = K)

    for (t in seq_along(x.b)) {
      tau <- get_tau(t, K, h, x.b, mu)
      upsilon <- get_upsilon(t, K, h, x.b, mu, kappa)

      # --- Conditional updates ---
      kappa.tmp[t, ] <- sqrt(kappa^2 + colSums(arcoef * tau)^2)
      mu.tmp[t, ] <- mu + atan(colSums(arcoef * upsilon))
    }

    kappa.t[[b]] <- kappa.tmp
    mu.t[[b]] <- mu.tmp
  }

  # --- Flatten lists to match original vector length ---
  kappa.t <- do.call(rbind, kappa.t)
  mu.t <- do.call(rbind, mu.t)

  # --- Q function matrix (weighted log-density) ---
  Q <- matrix(ncol = K, nrow = length(x))
  for (k in seq_len(K)) {
    for (i in seq_along(x)) {
      Q[i, k] <- w[i, k] *
        circular::dvonmises(x[i], mu.t[i, k], kappa.t[i, k], log = TRUE)
    }
  }

  -sum(Q) # return negative to maximize via optim
}

# --------------------------------------------------------------------------- #
#   MAIN FUNCTION (EM ALGORITHM)
# --------------------------------------------------------------------------- #
#' @description
#' fit circular mar model via em algorithm
#'
#' @param x circular data vector
#' @param burst vector identifying independent trajectories
#' @param K number of mixture components
#' @param h autoregressive lag order
#' @param tol convergence tolerance (delta log-likelihood)
#' @param maxit maximum em iterations
#' @return object of class cmar with estimated parameters
fit_cmar <- function(x, burst, K, h, tol = 1e-4, maxit = 100) {
  # init
  prob <- rep(1 / K, K)
  mu <- runif(K, 0, 2 * pi)
  lkappa <- log(runif(K, 1, 5))
  arcoef <- runif(h * K, 0, 0.5)

  param <- c(mu, lkappa, arcoef)

  l <- -Inf
  l[2] <- loglik(param, prob, x, burst, K, h)[[2]]
  it <- 2
  converged <- FALSE

  cat("\nStarting EM Algorithm for Circular MAR Model\n")

  # --- EM loop ---
  while (it < maxit && !converged) {
    # --- E-step ---
    f <- loglik(param, prob, x, burst, K, h)[[1]]
    w <- f / rowSums(f)

    # --- M-step ---
    prob <- colSums(w) / sum(colSums(w))

    init <- Sys.time()
    opt <- optim(
      par = param,
      fn = Q,
      x = x,
      w = w,
      burst = burst,
      K = K, h = h,
      control = list(maxit = 10000)
    )
    time <- Sys.time() - init

    if (opt$convergence != 0) { # convergence check
      warning("Iteration: ", it, "; Convergence: ", opt$convergence)
    }

    # --- Parameters update ...
    param <- opt$par
    it <- it + 1
    l[it] <- loglik(param, prob, x, burst, K, h)[[2]]

    cat(sprintf( # message
      "Iteration %d completed
      LogLik: %.4f
      Delta: %.6f
      Time M-step: %.2f sec\n
      ",
      it - 1,
      l[it],
      l[it] - l[it - 1],
      as.numeric(time, units = "secs")
    ))

    if (abs(l[it] - l[it - 1]) < tol) converged <- TRUE
  }

  # --- Post processing ---
  # transform log-kappa back to standard kappa
  param[(K + 1):(2 * K)] <- exp(param[(K + 1):(2 * K)]) # transform log-kappa

  names(param) <- c(
    paste0("mu", 1:K),
    paste0("kappa", 1:K),
    paste0("ar.l", rep(1:h, each = K), ".c", rep(1:K, times = h))
  )

  # --- Output object ---
  res <- list(
    coefficients = param,
    mixing_probs = prob,
    loglik       = l[it],
    history      = l,
    iterations   = it,
    convergence  = converged
  )

  class(res) <- "cmar"
  return(res)
}

# --------------------------------------------------------------------------- #
#   S3 METHODS (OUTPUT FORMATTING)
# --------------------------------------------------------------------------- #

#' print method for cmar objects
#'
#' @param x cmar object
#' @param ... further arguments passed to or from other methods
#' @export
print.cmar <- function(x, ...) {
  cat("\ncircular mixture autoregressive model\n")
  cat("-------------------------------------\n")

  cat("mixing probabilities:\n")
  m_probs <- matrix(x$mixing_probs, nrow = 1)
  colnames(m_probs) <- paste0("comp.", seq_along(x$mixing_probs))
  rownames(m_probs) <- "prob."
  print(round(m_probs, 4))

  cat("\ncoefficients:\n")
  # tabella pulita per la ui
  coef_df <- data.frame(
    parameter = names(x$coefficients),
    value     = round(x$coefficients, 4)
  )
  print(coef_df, row.names = FALSE, right = FALSE)

  cat("\nlog-likelihood:", round(x$loglik, 4), "\n")
}

#' summary method for cmar objects
#'
#' @param object cmar object
#' @param ... further arguments passed to or from other methods
#' @export
summary.cmar <- function(object, ...) {
  # param count: (K - 1) for weights + length of all other params
  K <- length(object$mixing_probs)
  n_param <- (K - 1) + length(object$coefficients)

  # aic calculation
  aic <- 2 * n_param - 2 * object$loglik

  # create summary object
  ans <- list(
    coefficients = object$coefficients,
    mixing_probs = object$mixing_probs,
    loglik       = object$loglik,
    iterations   = object$iterations,
    convergence  = object$convergence,
    aic          = aic,
    n_param      = n_param
  )

  class(ans) <- "summary.cmar"
  return(ans)
}

#' print method for summary.cmar objects
#'
#' @param x summary.cmar object
#' @param ... further arguments passed to or from other methods
#' @export
print.summary.cmar <- function(x, ...) {
  cat("\n==============================================\n")
  cat("   circular mixture autoregressive model\n")
  cat("==============================================\n")

  cat("\nmodel fit statistics:\n")
  cat(sprintf("log-likelihood: %.4f\n", x$loglik))
  cat(sprintf("aic:            %.4f\n", x$aic))
  cat(sprintf("parameters:     %d\n", x$n_param))

  cat("\nmixing probabilities:\n")
  m_probs <- matrix(x$mixing_probs, nrow = 1)
  colnames(m_probs) <- paste0("comp.", seq_along(x$mixing_probs))
  rownames(m_probs) <- "prob."
  print(round(m_probs, 4))

  cat("\nestimated coefficients:\n")
  # tabella pulita per la ui
  coef_df <- data.frame(
    parameter = names(x$coefficients),
    value     = round(x$coefficients, 4)
  )
  print(coef_df, row.names = FALSE, right = FALSE)

  cat("\nconvergence details:\n")
  cat("converged:  ", x$convergence, "\n")
  cat("iterations: ", x$iterations, "\n")
  cat("==============================================\n")
}
