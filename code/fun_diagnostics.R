pacman::p_load(
  circular,
  tidyverse,
  sf,
  patchwork,
  knitr
)
# --------------------------------------------------------------------------- #
#   CHARTS SIM
# --------------------------------------------------------------------------- #
circ.hist.sim <- function(dat) {
  ggplot(dat, aes(x = y, fill = factor(z))) +
    geom_histogram(
      breaks = seq(0, 2 * pi, length.out = 31),
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(0, 2 * pi),
      breaks = c(0, pi / 2, pi, 3 * pi / 2),
      labels = c("0", "π/2", "π", "3π/2")
    ) +
    coord_polar(
      theta = "x",
      direction = -1,
      start = 3 * pi / 2
    ) +
    theme_bw()
}

# Path
path.chart.sim <- function(dat, tburst, n) {
  dat <- dat %>%
    filter(burst == tburst) %>%
    mutate(
      D.x = cos(cumsum(x) %% (2 * pi)),
      D.y = sin(cumsum(x) %% (2 * pi)),
      xend = cumsum(D.x), yend = cumsum(D.y),
      x0 = lag(xend, default = 0),
      y0 = lag(yend, default = 0)
    ) %>%
    head(n)

  ggplot(dat, aes(
    x = x0, y = y0, xend = xend, yend = yend,
    color = as.factor(z)
  )) +
    geom_segment(
      arrow = arrow(length = unit(0.2, "cm")),
      linewidth = .75
    ) +
    coord_fixed() +
    labs(color = "Comp", title = "Real") +
    theme_bw()
}

# --------------------------------------------------------------------------- #
#   CHARTS FIT SIM
# --------------------------------------------------------------------------- #
# --- Functions ---
# Circular histogram
plot.polar.hist <- function(dat, y, z, title) {
  ggplot2::ggplot(dat, ggplot2::aes(x = y, fill = factor(z))) +
    ggplot2::geom_histogram(
      breaks = seq(0, 2 * pi, length.out = 31),
      color = "black"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 2 * pi),
      breaks = c(0, pi / 2, pi, 3 * pi / 2),
      labels = c("0", "π/2", "π", "3π/2")
    ) +
    coord_polar(
      theta = "x",
      direction = -1,
      start = 3 * pi / 2
    ) +
    labs(title = title, fill = "State", x = NULL, y = NULL) +
    theme_bw()
}

circ.hist <- function(sim, fit, mod) {
  # Label switching
  distmat <- outer(mod$mu, fit$params$mu, function(x, y) 1 - cos(x - y))
  idx <- integer()
  for (j in seq_len(mod$K)) {
    idx[j] <- which.min(distmat[j, ])
    distmat[, idx[j]] <- Inf
  }

  w <- fit$weights[, idx]

  dat <- sim %>%
    group_by(id, burst) %>%
    slice(-(1:mod$h)) %>%
    ungroup() %>%
    mutate(z.post = apply(w, 1, which.max))

  # parchwok
  p1 <- plot.polar.hist(dat, y, z, "True States")
  p2 <- plot.polar.hist(dat, y, z.post, "Estimated States")

  p <- p1 + p2 + patchwork::plot_layout(guides = "collect")

  return(p)
}

# Time series
time.series <- function(dat, tburst, n) {
  dat <- dat %>%
    filter(burst == tburst) %>%
    mutate(t = seq_along(x)) %>%
    head(n)

  p1 <- ggplot(dat, aes(x = t, y = x, color = factor(z))) +
    geom_point(size = 2) +
    geom_line(color = "black", linetype = "dashed") +
    scale_y_continuous(
      limits = c(0, 2 * pi),
      breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
      labels = c("0", "π/2", "π", "3π/2", "2π")
    ) +
    labs(color = "comp", title = "Real", y = "Angle") +
    theme_bw()

  p2 <- ggplot(dat, aes(x = t, y = x, color = factor(z.post))) +
    geom_point(size = 2) +
    geom_line(color = "black", linetype = "dashed") +
    scale_y_continuous(
      limits = c(0, 2 * pi),
      breaks = c(0, pi / 2, pi, 3 * pi / 2, 2 * pi),
      labels = c("0", "π/2", "π", "3π/2", "2π")
    ) +
    labs(color = "comp", title = "Estimated", y = "Angle") +
    theme_bw()

  p1 / p2
}

# Path
path.chart <- function(dat, tburst, n) {
  dat <- dat %>%
    filter(burst == tburst) %>%
    mutate(
      xend = cumsum(D.x), yend = cumsum(D.y),
      x0 = lag(xend, default = 0),
      y0 = lag(yend, default = 0)
    ) %>%
    head(n)

  p1 <- ggplot(dat, aes(
    x = x0, y = y0, xend = xend, yend = yend,
    color = as.factor(z)
  )) +
    geom_segment(
      arrow = arrow(length = unit(0.2, "cm")),
      linewidth = .75
    ) +
    coord_fixed() +
    labs(color = "Comp", title = "Real") +
    theme_bw()

  p2 <- ggplot(dat, aes(
    x = x0, y = y0, xend = xend, yend = yend,
    color = as.factor(z.post)
  )) +
    geom_segment(
      arrow = arrow(length = unit(0.2, "cm")),
      linewidth = .75
    ) +
    coord_fixed() +
    labs(color = "Comp", title = "Estimated") +
    theme_bw()

  p1 + p2
}

# Arrows time series
arrows.ts <- function(dat, tburst, n) {
  dat <- dat %>%
    filter(burst == tburst) %>%
    head(n) %>%
    mutate(
      t = seq_along(x),
      c = (t - 1) %% 5,
      r = (t - 1) %/% 5 + 1,
      x0 = 2 * c,
      y0 = 0,
      xend = x0 + D.x
    ) %>%
    pivot_longer(
      cols = c(z, z.post),
      names_to = "comp.type",
      values_to = "z"
    )

  ggplot(dat, aes(
    x = x0, y = y0, xend = xend, yend = D.y, color = factor(z)
  )) +
    geom_hline(yintercept = 0) +
    geom_segment(
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      linewidth = 0.8
    ) +
    coord_fixed(ratio = 1) +
    facet_grid(r ~ comp.type,
      labeller = labeller(comp.type = c("z" = "Real", "z.post" = "Est"))
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.y = element_blank(),
    ) +
    labs(color = "Comp")
}

# --------------------------------------------------------------------------- #
# SIMULATION STUDY ----
# --------------------------------------------------------------------------- #
# Model selection ----
mod.selection <- function(K, scenario) {
  load(paste0("data/simList_K", K, "n", scenario, ".RData"))
  n <- length(simlist)

  # Kest 2
  # WARNING: Sto tenendo a prescindere dalla convergenza
  load(paste0("data/fit_K", K, "Kest2n", scenario, ".RData"))
  convergedKest2 <- sapply(1:n, function(i) fitlist[[i]]$converged)
  bicKest2 <- sapply(1:n, function(i) fitlist[[i]]$criteria$bic)
  entropyKest2 <- sapply(1:n, function(i) fitlist[[i]]$criteria$entropy)
  iclKest2 <- bicKest2 + 2 * entropyKest2

  # Kest 3
  load(paste0("data/fit_K", K, "Kest3n", scenario, ".RData"))
  convergedKest3 <- sapply(1:n, function(i) fitlist[[i]]$converged)
  bicKest3 <- sapply(1:n, function(i) fitlist[[i]]$criteria$bic)
  entropyKest3 <- sapply(1:n, function(i) fitlist[[i]]$criteria$entropy)
  iclKest3 <- bicKest3 + 2 * entropyKest3

  # Kest 4
  load(paste0("data/fit_K", K, "Kest4n", scenario, ".RData"))
  convergedKest4 <- sapply(1:n, function(i) fitlist[[i]]$converged)
  bicKest4 <- sapply(1:n, function(i) fitlist[[i]]$criteria$bic)
  entropyKest4 <- sapply(1:n, function(i) fitlist[[i]]$criteria$entropy)
  iclKest4 <- bicKest4 + 2 * entropyKest4

  tibble(
    Kest = c(
      rep(2, n),
      rep(3, n),
      rep(4, n)
    ),
    converged = c(convergedKest2, convergedKest3, convergedKest4),
    bic = c(bicKest2, bicKest3, bicKest4),
    entropy = c(entropyKest2, entropyKest3, entropyKest4),
    icl = c(iclKest2, iclKest3, iclKest4)
  )

  tibble(
    replica = 1:n,
    converged = paste0(apply(
      cbind(convergedKest2, convergedKest3, convergedKest4), 1, sum
    ), "/3"),
    bic = apply(
      cbind(bicKest2, bicKest3, bicKest4), 1, which.min
    ) + 1,
    icl = apply(
      cbind(iclKest2, iclKest3, iclKest4), 1, which.min
    ) + 1
  )
}

# Parameter estimation ----
get.res.row <- function(i, fitlist, simlist, mod, formula) {
  fit <- fitlist[[i]]
  sim <- simlist[[i]] %>% mutate(burst = paste0("id", id, "burst", burst))

  # Label switching
  distmat <- outer(mod$mu, fit$params$mu, function(x, y) 1 - cos(x - y))
  idx <- integer()
  for (j in seq_len(mod$K)) { # get index order
    idx[j] <- which.min(distmat[j, ])
    distmat[, idx[j]] <- Inf
  }

  # order parameters
  mu <- fit$params$mu[idx]
  kappa <- fit$params$kappa[idx]
  arcoef <- fit$params$arcoef[, idx, drop = FALSE]

  # order weights
  w <- fit$weights[, idx]

  # order betas
  xmat <- model.matrix(formula, sim)
  xmat <- get.dat(sim$y, xmat, sim$burst, mod$h)$xmat
  formulamod <- as.formula("preds ~ . - 1")
  datmod <- data.frame(preds = I(w), xmat)
  fitGlm <- nnet::multinom(formulamod, data = datmod, trace = FALSE)
  betas <- coef(fitGlm)
  preds <- predict(fitGlm, type = "probs")

  # ARI
  sim <- sim %>%
    dplyr::arrange(burst, t) %>%
    dplyr::group_by(id, burst) %>%
    dplyr::slice(-(1:mod$h)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(z.post = apply(w, 1, which.max))

  # Error
  muErr <- (mu - mod$mu + pi) %% (2 * pi) - pi
  kappaErr <- kappa - mod$kappa
  arcoefErr <- arcoef - mod$arcoef
  betaErr <- betas - mod$beta.array

  # --- Dataset row ---
  resRow <- tibble(
    run = i,
    converged = fit$converged,
    ari = mclust::adjustedRandIndex(sim$z, sim$z.post)
  )

  # Circular and AR params
  for (k in seq_len(mod$K)) {
    resRow[[paste0("muEst_", k)]] <- mu[k]
    resRow[[paste0("muErr_", k)]] <- muErr[k]
    resRow[[paste0("kappaEst_", k)]] <- kappa[k]
    resRow[[paste0("kappaErr_", k)]] <- kappaErr[k]
    resRow[[paste0("arcoefEst_", k)]] <- arcoef[1, k] # Assumendo h=1
    resRow[[paste0("arcoefErr_", k)]] <- arcoefErr[1, k]
  }

  # Multinomial params (state 2:K)
  cov_names <- colnames(betas)
  for (k in 2:mod$K) {
    for (p in seq_along(cov_names)) {
      cname <- make.names(cov_names[p]) # Rimuove parentesi es. da (Intercept)
      resRow[[paste0("betaEst_St", k, "_", cname)]] <- betas[k - 1, p]
      resRow[[paste0("betaErr_St", k, "_", cname)]] <- betaErr[k - 1, p]
    }
  }

  resRow
}

res.summary <- function(res, mod) {
  res <- res %>% filter(converged == TRUE)
  ariMean <- mean(res$ari)

  # Init output list
  out <- list()

  # Circular and AR params
  for (k in 1:mod$K) {
    # Mu
    muEst <- res[[paste0("muEst_", k)]]
    muErr <- res[[paste0("muErr_", k)]]
    out[[length(out) + 1]] <- dplyr::tibble(
      Parameter = "mu", State = k, True = mod$mu[k],
      Est = as.numeric(
        suppressWarnings(circular::mean.circular(muEst))
      ),
      Bias = mean(muErr),
      SD = as.numeric(suppressWarnings(circular::angular.deviation(muEst))),
      RMSE = sqrt(mean(muErr^2))
    )

    # Kappa
    kappaEst <- res[[paste0("kappaEst_", k)]]
    kappaErr <- res[[paste0("kappaErr_", k)]]
    out[[length(out) + 1]] <- dplyr::tibble(
      Parameter = "kappa", State = k, True = mod$kappa[k],
      Est = mean(kappaEst),
      Bias = mean(kappaErr),
      SD = sd(kappaEst),
      RMSE = sqrt(mean(kappaErr^2))
    )

    # Arcoef
    arcoefEst <- res[[paste0("arcoefEst_", k)]]
    arcoefErr <- res[[paste0("arcoefErr_", k)]]
    out[[length(out) + 1]] <- dplyr::tibble(
      Parameter = "phi", State = k, True = mod$arcoef[1, k],
      Est = mean(arcoefEst),
      Bias = mean(arcoefErr),
      SD = sd(arcoefEst),
      RMSE = sqrt(mean(arcoefErr^2))
    )
  }

  # Betas
  cols <- grep("^betaEst_St2_", names(res), value = TRUE)
  cnames <- sub("^betaEst_St2_", "", cols)

  for (k in 2:mod$K) {
    for (p in seq_along(cnames)) {
      cname <- cnames[p]

      bEst <- res[[paste0("betaEst_St", k, "_", cname)]]
      bErr <- res[[paste0("betaErr_St", k, "_", cname)]]

      out[[length(out) + 1]] <- dplyr::tibble(
        Parameter = gsub("^X\\.|\\.$", "", cname),
        State = k,
        True = mod$beta.array[k - 1, p],
        Est = mean(bEst),
        Bias = mean(bErr),
        SD = sd(bEst),
        RMSE = sqrt(mean(bErr^2))
      )
    }
  }

  out <- dplyr::bind_rows(out)

  attr(out, "ARIMean") <- round(ariMean, 3)
  attr(out, "Converged") <- nrow(res)

  return(out)
}
