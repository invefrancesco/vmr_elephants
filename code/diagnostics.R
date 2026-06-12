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
circ.hist.sim(dat)
path.chart.sim(dat, 1, 50)

# --------------------------------------------------------------------------- #
#   CHARTS FIT
# --------------------------------------------------------------------------- #
# --- Functions ---
# Circular histogram
circ.hist <- function(dat) {
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

# --- Examples ---
load("data/sim6.RData")
load("data/fit6.RData")

# --- Data preparation ---
dat <- dat %>%
  group_split(burst) %>%
  map_dfr(\(df) df[-1, ]) %>%
  mutate(
    z.post = apply(fit$weights, 1, which.max),
    D.x = cos(x),
    D.y = sin(x)
  )

# ARI ----
mclust::adjustedRandIndex(dat$z, dat$z.post)

# Coef table ----
bind_rows(
  tibble(
    state = 1:mod$K,
    type = "Real",
    mu = mod$mu,
    kappa = mod$kappa,
    arcoef = as.numeric(mod$arcoef[1, ]), # (h = 1)
    prob = mod$prob
  ),
  tibble(
    state = seq_along(fit$params$mu),
    type = "Est",
    mu = fit$params$mu,
    kappa = fit$params$kappa,
    arcoef = as.numeric(fit$params$arcoef[1, ]),
    prob = fit$params$prob
  )
) %>%
  arrange(state, type) %>%
  knitr::kable(
    align = "c",
    digits = 3
  )

# Vusal ----
circ.hist(dat)
arrows.ts(datDgn5, 5, 50)
path.chart(datDgn5, 5, 50)
time.series(datDgn5, 5, 50)
