# Read command line arguments
args2 <- commandArgs(trailingOnly = TRUE)

# Load required packages
pacman::p_load(
  circular,
  tidyverse,
  doParallel,
  foreach
)

# Load the simulated data
name <- args2[1]
load(name)

# Extract parameters from arguments
K <- as.numeric(args2[2])
Kest <- as.numeric(args2[3])
n <- as.numeric(args2[4])
h <- 1

# Set up parallel backend
ncores <- 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Record start time
init <- Sys.time()

# Parallel loop over the number of simulations
fitlist <- foreach(
  i = 1:nSim,
  .export = c(".GlobalEnv"),
  .inorder = FALSE,
  .packages = c(
    "circular",
    "tidyverse",
    "sf"
  )
) %dopar% {
  # Log the start of the iteration to an external text file
  pid <- Sys.getpid()
  logFile <- paste0(".log/", pid, ".txt")
  sink(logFile, append = TRUE)

  cat(sprintf("Core %d: Starting fit for dataset %d of %d\n", Sys.getpid(), i, nSim),
    file = ".log/log.txt", append = TRUE
  )

  # Source estimation functions
  source("code/fun_estimation.R")

  # Prepare the dataset by creating a unique global burst ID
  dat <- simlist[[i]] %>% mutate(burst = paste0("id", id, "burst", burst))
  formula <- ~ env + sex

  # Fit the model
  set.seed(1234)
  fit <- cmar.ms(
    y = dat$y,
    burst = dat$burst,
    formula = formula,
    dat = dat,
    K = Kest, h = h,
    maxit = 1000
  )

  cat(sprintf("Core %d: Completed dataset %d!\n", Sys.getpid(), i),
    file = ".log/log.txt", append = TRUE
  )

  sink()
  return(fit)
}
stopCluster(cl)
print(Sys.time() - init)

# Save the final results
save(fitlist, file = paste0("data/fit_K", K, "Kest", Kest, "n", n, ".RData"))
