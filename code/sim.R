#' # circular mixture autoregressive (cmar) model
#'
#' this is the script for simulating data.
args <- commandArgs(trailingOnly = TRUE)

pacman::p_load(
  circular,
  tidyverse,
  sf,
  doParallel,
  foreach
)

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

nSim <- 20
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
seeds <- sample(1:10000, size = nSim, replace = FALSE)
ncores <- 3
cl <- makeCluster(ncores)
registerDoParallel(cl)
init <- Sys.time()
simlist <- foreach(
  i = 1:nSim,
  .export = c(".GlobalEnv"),
  .inorder = FALSE,
  .packages = c(
    "circular",
    "tidyverse",
    "sf"
  )
) %dopar% {
  source("code/fun_simulation.R")
  set.seed(seeds[i])
  sim.data(dat, mod, id, burst)
}
stopCluster(cl)
print(Sys.time() - init)

# Save ----
name <- paste0("data/simList_K", K, "n", args[length(args)], ".RData")
save.image(file = name)
