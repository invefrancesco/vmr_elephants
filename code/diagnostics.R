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
  ggplot(dat, aes(x = x, fill = factor(z))) +
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

# --- Examples ---
load("data/sim6.RData")
circ.hist.sim(simlist[[1]])
path.chart.sim(dat, 1, 50)

# --------------------------------------------------------------------------- #
#   CHARTS FIT SIM
# --------------------------------------------------------------------------- #
# --- Functions ---
# Circular histogram
circ.hist <- function(sim, fit) {
  dat <- sim %>%
    group_by(id, burst) %>%
    slice(-(1:mod$h)) %>%
    ungroup() %>%
    mutate(z.post = apply(fit$weights, 1, which.max))
  
  p1 <- ggplot(dat, aes(x = x, fill = factor(z))) +
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

  p2 <- ggplot(dat, aes(x = x, fill = factor(z.post))) +
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

  print(p1 + p2)
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
#   MODEL SELECTION
# --------------------------------------------------------------------------- #
load("data/simList_K2n1.RData")

# BIC ----
n <- length(fitlist)
("data/fit_K2Kest3n1.RData")
indx <- which(sapply(1:n, function(i) fitlist[[i]]$converged))

load("data/fit_K2Kest2n1.RData")
bicK2 <- sapply(indx, function(i) fitlist[[i]]$criteria$bic)

load("data/fit_K2Kest3n1.RData")
bicK3 <- sapply(indx, function(i) fitlist[[i]]$criteria$bic)

BIC <- apply(cbind(bicK2, bicK3), 1, which.min) + 1

# Entropy ----
load("data/fit_K2Kest2n1.RData")
eK2 <- sapply(
  indx,
  function(i) -sum(fitlist[[i]]$weights * log(fitlist[[1]]$weights + 1e-15))
)

load("data/fit_K2Kest3n1.RData")
eK3 <- sapply(
  indx,
  function(i) -sum(fitlist[[i]]$weights * log(fitlist[[1]]$weights + 1e-15))
)
# ICL ----
iclK2 <- bicK2 + 2 * eK2
iclK3 <- bicK3 + 2 * eK3
ICL <- apply(cbind(iclK2, iclK3), 1, which.min) + 1

sum(BIC == 2) / length(indx)
sum(ICL == 2) / length(indx)

# --------------------------------------------------------------------------- #
#   SIMULATION RESULTS
# --------------------------------------------------------------------------- #
get.res.row <- function(i, fitlist, simlist, mod) {
  fit <- fitlist[[i]]
  sim <- simlist[[i]]

  # Label switching
  distmat <- outer(mod$mu, fit$params$mu, function(x, y) 1 - cos(x - y))
  idx <- integer()
  for (j in seq_len(mod$K)) { # get index order
    idx[j] <- which.min(distmat[j, ])
    distmat[, idx[j]] <- Inf
  }

  # order parameters and weights
  mu <- fit$params$mu[idx]
  kappa <- fit$params$kappa[idx]
  arcoef <- fit$params$arcoef[, idx, drop = FALSE]
  prob <- fit$params$prob[idx]
  w <- fit$weights[, idx]

  # ARI
  sim <- sim %>%
    group_by(id, burst) %>%
    slice(-(1:mod$h)) %>%
    ungroup() %>%
    mutate(z.post = apply(w, 1, which.max))

  ari <- mclust::adjustedRandIndex(sim$z, sim$z.post)

  # Error
  muErr <- (mu - mod$mu + pi) %% (2 * pi) - pi
  kappaErr <- kappa - mod$kappa
  arcoefErr <- arcoef - mod$arcoef
  probErr <- prob - mod$prob

  # Dataset row
  resRow <- tibble(
    run = i,
    converged = TRUE,
    ari = ari
  )

  for (k in seq_len(mod$K)) {
    resRow[[paste0("muEst_", k)]] <- mu[k]
    resRow[[paste0("muErr_", k)]] <- muErr[k]

    resRow[[paste0("kappaEst_", k)]] <- kappa[k]
    resRow[[paste0("kappaErr_", k)]] <- kappaErr[k]

    resRow[[paste0("arcoefEst_", k)]] <- arcoef[1, k]
    resRow[[paste0("arcoefErr_", k)]] <- arcoefErr[k]

    resRow[[paste0("probEst_", k)]] <- prob[k]
    resRow[[paste0("probErr_", k)]] <- probErr[k]
  }

  resRow
}

res.summary <- function(res, mod) {
  ariMean <- mean(res$ari, na.rm = TRUE)
  ariSD <- sd(res$ari, na.rm = TRUE)

  muMean <- numeric(mod$K)
  muSD <- numeric(mod$K)
  muBias <- numeric(mod$K)
  muRMSE <- numeric(mod$K)

  kappaMean <- numeric(mod$K)
  kappaSD <- numeric(mod$K)
  kappaBias <- numeric(mod$K)
  kappaRMSE <- numeric(mod$K)

  arcoefMean <- numeric(mod$K)
  arcoefSD <- numeric(mod$K)
  arcoefBias <- numeric(mod$K)
  arcoefRMSE <- numeric(mod$K)

  probMean <- numeric(mod$K)
  probSD <- numeric(mod$K)
  probBias <- numeric(mod$K)
  probRMSE <- numeric(mod$K)

  for (k in 1:mod$K) {
    # --- mu ---
    colEst <- paste0("muEst_", k)
    colErr <- paste0("muErr_", k)

    suppressWarnings({
      muMean[k] <- as.numeric(circular::mean.circular(res[[colEst]]))
      muSD[k] <- as.numeric(circular::angular.deviation(res[[colEst]]))
    })
    muBias[k] <- mean(res[[colErr]], na.rm = TRUE)
    muRMSE[k] <- sqrt(mean(res[[colErr]]^2, na.rm = TRUE))

    # --- kappa ---
    colEst <- paste0("kappaEst_", k)
    colErr <- paste0("kappaErr_", k)

    kappaMean[k] <- mean(res[[colEst]], na.rm = TRUE)
    kappaSD[k] <- sd(res[[colEst]], na.rm = TRUE)
    kappaBias[k] <- mean(res[[colErr]], na.rm = TRUE)
    kappaRMSE[k] <- sqrt(mean(res[[colErr]]^2, na.rm = TRUE))

    # --- arcoef ---
    colEst <- paste0("arcoefEst_", k)
    colErr <- paste0("arcoefErr_", k)

    arcoefMean[k] <- mean(res[[colEst]], na.rm = TRUE)
    arcoefSD[k] <- sd(res[[colEst]], na.rm = TRUE)
    arcoefBias[k] <- mean(res[[colErr]], na.rm = TRUE)
    arcoefRMSE[k] <- sqrt(mean(res[[colErr]]^2, na.rm = TRUE))

    # --- prob ---
    colEst <- paste0("probEst_", k)
    colErr <- paste0("probErr_", k)

    probMean[k] <- mean(res[[colEst]], na.rm = TRUE)
    probSD[k] <- sd(res[[colEst]], na.rm = TRUE)
    probBias[k] <- mean(res[[colErr]], na.rm = TRUE)
    probRMSE[k] <- sqrt(mean(res[[colErr]]^2, na.rm = TRUE))
  }

  out <- tibble::tibble(
    state = 1:mod$K,
    muMean = muMean, muSD = muSD, muBias = muBias, muRMSE = muRMSE,
    kappaMean = kappaMean, kappaSD = kappaSD, kappaBias = kappaBias, kappaRMSE = kappaRMSE,
    arcoefMean = arcoefMean, arcoefSD = arcoefSD, arcoefBias = arcoefBias, arcoefRMSE = arcoefRMSE,
    probMean = probMean, probBias = probBias, probRMSE = probRMSE
  )

  attr(out, "ARI_Mean") <- ariMean
  attr(out, "ARI_SD") <- ariSD

  return(out)
}

load("data/simList_K2n1.RData")
load("data/fit_K2Kest2n1.RData")

res <- map_dfr(1:20, function(i) {
  get.res.row(i = i, fitlist = fitlist, simlist = simlist, mod = mod)
})

res.summary(res, mod)
select(res, muEst_1, muEst_2)
