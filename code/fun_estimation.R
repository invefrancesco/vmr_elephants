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
# Function to get lags and target vector y.t
get.dat <- function(y, xmat, burst, h) {
  lag <- list()
  y.t <- list()
  xmat.t <- list()
  for (b in unique(burst)) {
    y.b <- y[burst == b]
    xmat.b <- xmat[burst == b, ]
    foo <- sapply(1:h, function(g) head(c(rep(NA, g), y.b), n = length(y.b)))

    lag[[as.character(b)]] <- foo[complete.cases(foo), , drop = FALSE]
    y.t[[as.character(b)]] <- y.b[complete.cases(foo)]
    xmat.t[[as.character(b)]] <- xmat.b[complete.cases(foo), ]
  }
  list(
    lag = do.call(rbind, lag),
    y.t = unlist(y.t, use.names = FALSE),
    xmat = do.call(rbind, xmat.t)
  )
}

# Function to transform natural to working parameters
n2w <- function(mu, kappa, arcoef) c(mu, log(kappa), arcoef)

# Function to transform working to natural parameters
w2n <- function(wpar, K, h) {
  mu <- wpar[1:K] %% (2 * pi)
  kappa <- exp(wpar[(K + 1):(2 * K)])
  arcoef <- matrix(tanh(wpar[(2 * K + 1):(2 * K + (h * K))]), nrow = h, ncol = K)
  list(mu = mu, kappa = kappa, arcoef = arcoef)
}

# --------------------------------------------------------------------------- #
#   LIKELIHOOD FUNCTION
# --------------------------------------------------------------------------- #
llk <- function(y.t, lag, zzz) {
  K <- length(zzz$mu)

  # time dependent mu and k
  kappa.t <- matrix(ncol = K, nrow = nrow(lag))
  mu.t <- matrix(ncol = K, nrow = nrow(lag))
  for (k in 1:K) {
    foo <- sin(lag - zzz$mu[k]) %*% zzz$arcoef[, k]
    kappa.t[, k] <- sqrt(zzz$kappa[k]^2 + foo^2)
    mu.t[, k] <- zzz$mu[k] + atan(foo / zzz$kappa[k])
  }
  dens <- (exp(kappa.t * (cos(y.t - mu.t) - 1)) /
    (2 * pi * besselI(kappa.t, 0, expon.scaled = TRUE)))

  # -llk
  list(
    dens = zzz$preds * dens,
    llk = sum(log(rowSums(zzz$preds * dens)))
  )
}

# --------------------------------------------------------------------------- #
#   Q FUNCTION
# --------------------------------------------------------------------------- #
Q <- function(wpar, y.t, lag, w, K, h) {
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
  ldens <- kappa.t * (cos(y.t - mu.t) - 1) -
    log(2 * pi * besselI(kappa.t, 0, expon.scaled = TRUE))

  -sum(w * ldens)
}
# --------------------------------------------------------------------------- #
#   MAIN FUNCTION
# --------------------------------------------------------------------------- #

cmar <- function(y, burst, formula, dat, K, h, tol = 1e-4, maxit = 500, verbose = TRUE) {
  # Init data
  xmat <- model.matrix(formula, dat)
  dat <- get.dat(y, xmat, burst, h)
  y.t <- dat$y.t
  lag <- dat$lag
  xmat <- dat$xmat

  # Init alpha and beta
  preds <- LaplacesDemon::rdirichlet(length(y.t), alpha = rep(1 / K, K))
  betas <- NULL

  # Compute first step alpha and beta
  datmod <- data.frame(preds = I(preds), xmat)
  formulamod <- as.formula("preds ~ . - 1")
  fit <- nnet::multinom(formulamod, data = datmod, trace = FALSE)
  betas <- coef(fit)
  preds <- predict(fit, type = "probs")

  # Init parameters
  zzz <- list(
    mu = runif(K, 0, 2 * pi),
    kappa = runif(K, 0, 10),
    arcoef = matrix(runif(h * K, 0, .5), nrow = h, ncol = K),
    betas = betas,
    preds = preds
  )

  # Init likelihood trace
  l <- -Inf
  l[2] <- llk(y.t, lag, zzz)$llk
  it <- 2
  converged <- FALSE

  # --- EM loop ---
  while (it < maxit && !converged) {
    init <- Sys.time()
    # --- E-step ---
    dens <- llk(y.t, lag, zzz)$dens
    w <- dens / rowSums(dens)

    # --- M-step ---
    # Betas update
    datmod <- data.frame(preds = I(w), xmat)
    fit <- nnet::multinom(formulamod, data = datmod, trace = FALSE)
    betas <- coef(fit)
    preds <- predict(fit, type = "probs")

    # Thetas update
    wpar <- n2w(zzz$mu, zzz$kappa, zzz$arcoef)
    opt <- optim(
      par = wpar,
      fn = Q,
      y.t = y.t,
      lag = lag,
      w = w,
      K = K, h = h,
      method = "BFGS",
      control = list(maxit = 10000)
    )

    # --- Parameters update ---
    zzz <- w2n(opt$par, K, h)
    zzz$betas <- betas
    zzz$preds <- preds
    it <- it + 1
    l[it] <- llk(y.t, lag, zzz)$llk

    if (verbose) {
      cat(sprintf( # Message
        "\nIteration %d completed | Delta llk: %.6f | Time: %.2f sec",
        it - 1,
        l[it] - l[it - 1],
        as.numeric(Sys.time() - init, units = "secs")
      ))
    }

    if (abs(l[it] - l[it - 1]) < tol) converged <- TRUE
  }

  # --- Output ---
  nobs <- length(y.t)
  npar <- K + K + (h * K) + length(zzz$betas)
  lhat <- l[it]

  obj <- list(
    params = zzz, # Parameters
    llk = lhat, # Log-likelihood
    iter = it - 1, # # of iterations
    converged = converged, # True if converged
    trace = l[2:it], # Log-likelihood trace
    weights = w, # Final weights
    criteria = list(
      aic = -2 * lhat + 2 * npar,
      bic = -2 * lhat + log(nobs) * npar,
      npar = npar,
      entropy = -sum(w * log(w + 1e-15))
    ),
    call = list(K = K, h = h, tol = tol)
  )

  class(obj) <- "cmar"
  obj
}

cmar.ms <- function(y, burst, formula, dat, K, h, tol = 1e-4, maxit = 500, starts = 5, verbose = TRUE) {
  init.ms <- Sys.time()
  fit <- list()
  for (i in 1:starts) {
    if (verbose) cat(sprintf("\n START: %d/%d \n", i, starts))

    foo <- try(
      cmar(
        y, burst,
        formula = formula, dat = dat, K = K, h = h,
        tol = tol, maxit = maxit, verbose = verbose
      ),
      silent = TRUE
    )

    if (inherits(foo, "try-error")) {
      fit[[i]] <- NULL
    } else {
      fit[[i]] <- foo
    }
  }
  if (verbose) {
    cat(sprintf(
      "\nTotal time: %.2f min",
      as.numeric(Sys.time() - init.ms, units = "mins")
    ))
  }
  fit <- fit[!sapply(fit, is.null)]
  fit[[which.max(sapply(fit, function(f) f$llk))]]
}
