#' # circular mixture autoregressive (cmar) model
#'
#' this is the main script.

pacman::p_load(
  circular,
  tidyverse,
  sf
)
source("code/fun_estimation.R")
source("code/fun_simulation.R")

# --------------------------------------------------------------------------- #
#   TRIAL PARAMETERS AND SIMULATION
# --------------------------------------------------------------------------- #

# --- design from real data ---
load("data/elephants.RData")

dat <- data %>%
  arrange(id, burst_, t2_) %>%
  mutate(id = consecutive_id(id)) %>%
  mutate(burst = consecutive_id(burst_), .by = id) %>%
  select(id, burst)

# --- Trial parameters ---
K <- 2
h <- 1

mod1 <- list(
  K = 2, h = 1,
  mu = c(0, pi),
  kappa = c(5, 5),
  arcoef = matrix(c(0.5, -0.5), nrow = 1, ncol = 2),
  prob = c(0.5, 0.5)
)

mod2 <- list(
  K = 2, h = 1,
  mu = c(0, pi / 4),
  kappa = c(1, 1.5),
  arcoef = matrix(c(0.2, 0.2), nrow = 1, ncol = 2),
  prob = c(0.7, 0.3)
)

mod3 <- list(
  K = 2, h = 1,
  mu = c(pi / 2, 3 * pi / 2),
  kappa = c(4, 4),
  arcoef = matrix(c(0.1, -0.8), nrow = 1, ncol = 2),
  prob = c(0.4, 0.6)
)

mod4 <- list(
  K = 2, h = 1,
  mu = c(0, pi / 2),
  kappa = c(2, 3),
  arcoef = matrix(c(0.1, 0.7), nrow = 1, ncol = 2),
  prob = c(0.3, 0.7)
)

mod5 <- list(
  K = 2, h = 1,
  mu = c(0, pi / 2),
  kappa = c(2, 3),
  arcoef = matrix(c(0.1, 0.7), nrow = 1, ncol = 2),
  prob = c(0.5, 0.5)
)
# --- sim data ---
set.seed(1234)
dat1 <- sim.data(dat, mod1, id, burst)
dat2 <- sim.data(dat, mod2, id, burst)
dat3 <- sim.data(dat, mod3, id, burst)
dat4 <- sim.data(dat, mod4, id, burst)
dat5 <- sim.data(dat, mod5, id, burst)

# --------------------------------------------------------------------------- #
#   ESTIMATION
# --------------------------------------------------------------------------- #
fit1 <- cmar.ms(dat1$x, dat1$burst, K, h)
fit2 <- cmar.ms(dat2$x, dat2$burst, K, h) # Did not converge
fit3 <- cmar.ms(dat3$x, dat3$burst, K, h)
fit4 <- cmar.ms(dat4$x, dat4$burst, K, h)
fit5 <- cmar.ms(dat5$x, dat5$burst, K, h)
