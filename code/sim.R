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

set.seed(1234)
dat <- tibble(
  id = rep(1:5, each = 1000),
  burst = rep(1:25, each = 200),
  env = rnorm(5000, 0, 1)
) %>% 
  group_by(id) %>% 
  mutate(sex = sample(0:1, 1, replace = TRUE)) %>% 
  ungroup()

# Trial parameters ----
parse_num <- function(x) {
  unname(sapply(x, function(val) eval(parse(text = val))))
}

nSim <- 4
h <- 1
K <- parse_num(args[1])
mu <- parse_num(args[2:(K + 1)])
kappa <- parse_num(args[(K + 2):(2 * K + 1)])
arcoef <- matrix(
  parse_num(args[(2 * K + 2):(2 * K + h * K + 1)]),
  nrow = h, ncol = K
)

if (K == 2) {
  beta.array <- matrix(
    c(
      .2, # Intercept: State 1 (Explore) is less likely
      -.8, # Env: Higher env more explore
      -.2 # Sex
    ),
    nrow = 1, ncol = 3, byrow = TRUE
  )
}
if (K == 3) {
  beta.array <- matrix(c(
    # Int, Env, Sex
    .2, .5, -.05, # 2 vs 1: positive env, negative sex
    .3, -.2, .8 # 3 vs 1: negative env, positive sex
  ), nrow = 2, ncol = 3, byrow = TRUE)
}
if (K == 4) {
  beta.array <- matrix(c(
    # Int, Env, Sex
    .5, .5, .01, # State 2
    .1, .1, -.8, # State 3
    -1, -.8, 0.5 # State 3
  ), nrow = 3, ncol = 3, byrow = TRUE)
}

# Define the model ----
mod <- list(
  K = K, h = h,
  mu = mu, kappa = kappa,
  arcoef = arcoef,
  beta.array = beta.array
)
print(mod)

# Simulate data ----
seeds <- sample(1:10000, size = nSim, replace = FALSE)
ncores <- 4
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
  sim.data(
    dat, mod,
    formula = ~ env + sex,
    id, burst
  )
}
stopCluster(cl)
print(Sys.time() - init)

# Save ----
name <- paste0("data/simList_K", K, "n", args[length(args)], ".RData")
save.image(file = name)
