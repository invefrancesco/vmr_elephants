#' # circular mixture autoregressive (cmar) model
#'
#' this is the script for simulating data.
args <- commandArgs(trailingOnly = TRUE)

pacman::p_load(
  circular,
  tidyverse,
  sf
)
# source("code/fun_estimation.R")
source("code/fun_simulation.R")

# --------------------------------------------------------------------------- #
#   TRIAL PARAMETERS AND SIMULATION
# --------------------------------------------------------------------------- #

# Design from real data ----
load("data/elephants.RData")

dat <- data %>%
  arrange(id, burst_, t2_) %>%
  mutate(id = consecutive_id(id)) %>%
  mutate(burst = consecutive_id(burst_), .by = id) %>%
  select(id, burst)

# Trial parameters ----
parse_num <- function(x) {
  unname(sapply(x, function(val) eval(parse(text = val))))
}

h <- 1
K <- parse_num(args[1])
mu <- parse_num(args[2:(K + 1)])
kappa <- parse_num(args[(K + 2):(2 * K + 1)])
arcoef <- matrix(
  parse_num(args[(2 * K + 2):(2 * K + h * K + 1)]),
  nrow = h, ncol = K
)
prob <- parse_num(args[(2 * K + h * K + 2):(3 * K + h * K + 1)])
prob <- prob / sum(prob)

# Define the model ----
mod <- list(K = K, h = h, mu = mu, kappa = kappa, arcoef = arcoef, prob = prob)
print(mod)

# Simulate data ----
set.seed(1234 + parse_num(args[length(args)]))
dat <- sim.data(dat, mod, id, burst)
summary(dat)

# Estimation ----
# fit <- cmar.ms(dat$x, dat$burst, K, h)

# Save ----
name <- paste0("data/sim", args[length(args)], ".RData")
save(mod, dat, file = name)
