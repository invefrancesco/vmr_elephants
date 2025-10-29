#' This file contains the code and the functions needed to estimate a mixture
#' of Von Mises on elephant data. 
#' 
#' The setting is the following: 
#' 
#' - We observe $Y_1, Y_2, \dots, Y_n$, a random sample where 
#' $Y_i \sim \text{Mixture}(\mu, \kappa, \pi)$.
#' - The mixture components distibutions are Von Mises
#' - The $\mu_k$ for each component are given by $2 \cdot \atan (\Beta X)$
#' 
#' ----------------------------------------------------------------------------
#' # Packages and data load
rm(list = ls())
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse, # Load multiple 'tidyverse' packages in a single step
  sf, # Support for simple feature access
  amt, # Manage and analyze animal movement data
  lubridate, # Functions to work with date-times and time-spans
  knitr, # Provides a general-purpose tool for dynamic report generation
  kableExtra # Build complex HTML or 'LaTeX' tables
)

dir_data <- paste0(getwd(), "/data")
load(paste0(dir_data, "/elephants.RData"))

#' Design matrix 
formula <- ~ distriv_std + elev_std + ndvi_std + seas + sex

data <- data %>%
  dplyr::select(all.vars(formula), "ta_", "burst_", "id") %>%
  tidyr::drop_na()

X <-  model.matrix(formula, data)

y <- data$ta_
# -----------------------------------------------------------------------------
#' # Functions
#' 
#' Function to be maximized in the M-step
m <- function(w, y, X, par){
  # Parameters
  # A matrix of betas k*B with k the number of mixture components and B the
  # number of covariates
  beta <- matrix(par[1:ncol(X)*ncol(w)], nrow = ncol(X), ncol = ncol(w))
  # A vector of k
  kappa <- par[(ncol(X)*ncol(w) + 1): (ncol(X)*ncol(w) + ncol(w))]
  # A matrix of observation-related means
  mu <- 2 * atan(X %*% beta)
  
  # Compute the value of the complete-data log-likelihood
  Q <- sapply(
    seq(ncol(w)),
    function(k)
      w[,k] * (kappa[k] * cos(y - mu[,k]) - log(2*pi*besselI(kappa[k], 0)))
  )
  -sum(Q)
}

# Log-likelihod function to store compute wieights and store the value
loglik <- function(p, y, X, par){
  # Parameters
  beta <- matrix(par[1:ncol(X)*ncol(w)], nrow = ncol(X), ncol = ncol(w))
  kappa <- par[(ncol(X)*ncol(w) + 1): (ncol(X)*ncol(w) + ncol(w))]
  mu <- 2 * atan(X %*% beta)
  
  # Compute the value of each component of the log-likelihood
  # It must be a matrix n * k
  l <- sapply(
    seq_along(p),
    function(k) 
      p[k] * (exp(kappa[k] * cos(y - mu[k])) / (2 * pi * besselI(kappa[k], 0)))
  )
  list(l, sum(log(rowSums(l))))
}

# EM-algorith,
  # E step
  f <- loglik(p, y, X, par)[[1]]
  w <- f/rowSums(f)
  
  # M-step 
  p <- sapply(seq_along(p), function(k) mean(w[, k]))
  
  par <- par
  opt <- optim(fn = m, w, y, X, par)
  
  # Update iteration value
  it <- it +1
  l[it] <- loglik(p, y, X, opt$par)