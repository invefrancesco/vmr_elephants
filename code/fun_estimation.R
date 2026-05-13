#' # circular mixture autoregressive (cmar) model
#'
#' this script implements the expectation-maximization (em) algorithm
#' to estimate parameters of a circular mar model using von mises distributions.
#' it handles independent bursts of circular time-series data.

pacman::p_load(
  circular
)

# --------------------------------------------------------------------------- #
#   AUXILIARY FUNCTIONS
# --------------------------------------------------------------------------- #
# Function to get lags and target vector x.t
get.dat <- function(x, burst, h) {
  lag <- list()
  x.t <- list()
  for (b in unique(burst)) {
    x.b <- x[burst == b]
    foo <- sapply(1:h, function(g) head(c(rep(NA, g), x.b), n = length(x.b)))

    lag[[as.character(b)]] <- foo[complete.cases(foo), , drop = FALSE]
    x.t[[as.character(b)]] <- x.b[complete.cases(foo)]
  }
  list(
    lag = do.call(rbind, lag),
    x.t = unlist(x.t)
  )
}

# Function to transform natural to working parameters
n2w <- function(mu, kappa, arcoef) c(mu, log(kappa), arcoef)

# Function to transform working to natural parameters
w2n <- function(wpar, K, h) {
  mu <- wpar[1:K] %% (2 * pi)
  kappa <- exp(wpar[(K + 1):(2 * K)])
  arcoef <- matrix(wpar[(2 * K + 1):(2 * K + (h * K))], nrow = h, ncol = K)
  list(mu = mu, kappa = kappa, arcoef = arcoef)
}
# --------------------------------------------------------------------------- #
#   LIKELIHOOD FUNCTION
# --------------------------------------------------------------------------- #
llk <- function(x.t, lag, zzz) {
  K <- length(zzz$mu)

  # time dependent mu and k
  kappa.t <- matrix(ncol = K, nrow = nrow(lag))
  mu.t <- matrix(ncol = K, nrow = nrow(lag))
  for (k in 1:K) {
    foo <- sin(lag - zzz$mu[k]) %*% zzz$arcoef[, k]
    kappa.t[, k] <- sqrt(zzz$kappa[k]^2 + foo^2)
    mu.t[, k] <- zzz$mu[k] + atan(foo / zzz$kappa[k])
  }
  dens <- (exp(kappa.t * (cos(x.t - mu.t) - 1)) /
    (2 * pi * besselI(kappa.t, 0, expon.scaled = TRUE)))

  # -llk
  list(
    dens = dens,
    llk = sum(log(dens %*% zzz$prob))
  )
}

# --------------------------------------------------------------------------- #
#   Q FUNCTION
# --------------------------------------------------------------------------- #
Q <- function(x.t, lag, wpar, w, K, h) {
  # --- Unpack parameters ---
  zzz <- w2n(wpar, K, h)

  # time dependent mu and k
  kappa.t <- matrix(ncol = K, nrow = nrow(lag))
  mu.t <- matrix(ncol = K, nrow = nrow(lag))
  for (k in 1:K) {
    foo <- sin(lag - zzz$mu[k]) %*% zzz$arcoef[, k]
    kappa.t[, k] <- sqrt(zzz$kappa[k]^2 + foo^2)
    mu.t[, k] <- zzz$mu[k] + atan(foo / zzz$kappa[k])
  }

  # --- Q function matrix (weighted log-density) ---
  # log(f) = kappa * (cos(x - mu) - 1) - log(2 * pi * besselI)
  ldens <- kappa.t * (cos(x.t - mu.t) - 1) -
    log(2 * pi * besselI(kappa.t, 0, expon.scaled = TRUE))

  -sum(w * ldens)
}
# --------------------------------------------------------------------------- #
#   MAIN FUNCTION
# --------------------------------------------------------------------------- #

cmar <- function(x, burst, K, h, tol = 1e-4, maxit = 100) {
  # Init data
  dat <- get.dat(x, burst, h)
  x.t <- dat$x.t
  lag <- dat$lag

  # Init parameters
  zzz <- list(
    mu = runif(K, 0, 2 * pi),
    kappa = runif(K, 0, 10),
    arcoef = matrix(runif(h * K, 0, .5), nrow = h, ncol = K),
    prob = rep(1 / K, K)
  )

  l <- -Inf
  l[2] <- llk(x.t, lag, zzz)$llk
  it <- 2
  converged <- FALSE

  cat("\nStarting EM Algorithm for Circular MAR Model\n")

  # --- EM loop ---
  while (it < maxit && !converged) {
    # --- E-step ---
    dens <- llk(x.t, lag, zzz)$dens
    w <- (dens * rep(zzz$prob, each = nrow(dens))) /
      rowSums(dens * rep(zzz$prob, each = nrow(dens)))

    # --- M-step ---
    wpar <- n2w(zzz$mu, zzz$kappa, zzz$arcoef)
    init <- Sys.time()
    opt <- optim(
      par = wpar,
      fn = Q,
      x.t = x.t,
      lag = lag,
      w = w,
      K = K, h = h,
      method = "BFGS",
      control = list(maxit = 10000)
    )
    time <- Sys.time() - init

    if (opt$convergence != 0) { # Convergence check
      warning("Iteration: ", it, "; Convergence: ", opt$convergence)
    }

    # --- Parameters update ---
    zzz <- w2n(opt$par, K, h)
    zzz$prob <- colSums(w) / sum(colSums(w))
    it <- it + 1
    l[it] <- llk(x.t, lag, zzz)$llk

    cat(sprintf( # Message
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

  # --- Output ---
  n.obs <- length(x.t)
  # K (mu) + K (kappa) + h*K (ar) + (K-1) (prob)
  n.par <- K + K + (h * K) + (K - 1)
  l.hat <- l[it]

  obj <- list(
    params = zzz, # Parametri naturali finali
    llk = l.hat, # Log-likelihood finale
    iter = it - 1, # Numero di iterazioni effettuate
    trace = l[2:it], # Storia della verosimiglianza per check convergenza
    weights = w, # Responsabilità (per clustering/segmentazione)
    criteria = list(
      aic = -2 * l.hat + 2 * n.par,
      bic = -2 * l.hat + log(n.obs) * n.par,
      npar = n.par
    ),
    call = list(K = K, h = h, tol = tol)
  )

  class(obj) <- "cmar"
  obj
}
