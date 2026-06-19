args2 <- commandArgs(trailingOnly = TRUE)

pacman::p_load(
  circular,
  tidyverse,
  sf,
  doParallel,
  foreach
)

args = c("data/simList.RData", 2, 2)

name <- args2[1]
load(name)

K <- as.numeric(args2[2])
Kest <- as.numeric(args2[3])
h <- 1

ncores <- 3
cl <- makeCluster(ncores)
registerDoParallel(cl)
init <- Sys.time()
fitlist <- foreach(
  i = 1:nSim,
  .export = c(".GlobalEnv"),
  .inorder = FALSE,
  .packages = c("circular",
                "tidyverse",
                "sf")
)%dopar%{
  source("code/fun_estimation.R")
  dat <- simlist[[i]] %>% mutate(burst = paste0("id", id, "burst", burst))
  set.seed(1234)
  fit <- cmar.ms(dat$x, dat$burst, Kest, h)
  return(fit)
}
stopCluster(cl)
print(Sys.time() - init)


# save(fitlist, file = paste0("data/fit", name, ".RData")) # Decidere nome

###############################################################################
# REAL DATA
###############################################################################
# 
# load("data/elephants.RData")
# 
# dat <- data %>%
#   mutate(
#     burst = paste0("id", id, "burst", burst_),
#     x = ta_ %% (2 * pi)
#   ) %>%
#   group_by(id, burst_) %>%
#   slice(-1) %>%
#   ungroup()
# 
# h <- 1
# fitK2 <- cmar.ms(dat$ta_, dat$burst, 2, h)
# fitK3 <- cmar.ms(dat$ta_, dat$burst, 3, h)
# fitK4 <- cmar.ms(dat$ta_, dat$burst, 4, h)
# save(fitK2, file = "data/fitK2.RData")

# TODO
# [ ] 
