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
#' # PACKAGES
#' ----------------------------------------------------------------------------
rm(list = ls())
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse, # Load multiple 'tidyverse' packages in a single step
  sf, # Support for simple feature access
  amt, # Manage and analyze animal movement data
  lubridate, # Functions to work with date-times and time-spans
  knitr, # Provides a general-purpose tool for dynamic report generation
  kableExtra, # Build complex HTML or 'LaTeX' tables
  collapse
)

#' ----------------------------------------------------------------------------
#' # FUNCTIONS
#' ----------------------------------------------------------------------------
#' - Auxiliary function to compute the lagged residuals by taking into account
#' id and burst.
get_lagged_res <- function(y, eta, burst, id) {
  K <- ncol(eta)
  y_matrix <- matrix(y, nrow = length(y), ncol = K)

  r <- as_tibble(y_matrix - eta) %>%
    fmutate(id = id, 
            burst = burst) %>% 
    fgroup_by(id, burst) %>% 
    fmutate(across(-c(id, burst), lag)) %>% 
    ungroup() %>% 
    select(-c(id, burst)) %>%
    qM()
}
#' ---------------------------------------------------------------------------- 
#' - Auxiliary function to compute a column of Q values
get_q <- function(k, w, y, mu_t, kappa_t) {
  fk <- (kappa_t[,k] * cos(y - mu_t[,k]) - log(2*pi*besselI(kappa_t[,k], 0)))
  w[,k] * fk
}
#' ----------------------------------------------------------------------------
#' - Function `m` to be maximized in the M-step
#' @param w A matrix n*K of observation-specific posterior weights.
#' @param y The response vector.
#' @param X The design matrix
#' @param par The vector of parameters to be maximized.
#' @param burst The vector of bursts.
#' @param id The factor of individual id.
m <- function(w, y, X, par, burst, id) {
  # Parameters
  K <- ncol(w) # The number of mixture components
  B <- ncol(X) # The number of covariates
  beta <- matrix(par[1:(B * K)], nrow = B, ncol = K) # B*K betas
  phi <- par[((B * K + 1) : (B * K + K))] # K auto-regressive parameters
  kappa <- exp(par[(B * K + K + 1) : (B * K + K + K)]) # K concentrations

  # Function
  eta <- 2 * atan(X %*% beta) %>% unname() # linear predictor
  r <- get_lagged_res(y, eta, burst, id) # lagged residuals
  kappa_t <- sqrt(kappa^2 + (phi * sin(r))^2) # dynamic kappa
  mu_t <- eta + atan( (phi * sin (r)) / kappa_t) # dynamic mu
  
  Q <- sapply(seq_len(K),
         get_q,
         w = w, 
         y = y,
         mu_t = mu_t,
         kappa_t = kappa_t)
  
  if(nrow(r[!complete.cases(r),]) != nrow(Q[!complete.cases(Q),])){
    warning(
      "# of groups:", nrow(r[!complete.cases(r),]), ", ",
      "# of NAs:", nrow(Q[!complete.cases(Q),]))
  }

  -sum(Q, na.rm = TRUE)
}
#' ----------------------------------------------------------------------------
#' - Auxiliary function to compute a column of likelihood values
get_L <- function(k, p, y, mu_t, kappa_t) {
  numerator <- exp(kappa_t[,k] * cos(y - mu_t[,k]))
  denominator <- 2 * pi * besselI(kappa_t[,k], 0)
  fk <- numerator / denominator
  
  p[k] * fk
}
#' - Log-likelihood function to compute weights and store the value
#' @param p a vector of mixture component's probabilities
#' @param y response vector.
#' @param X design matrix.
#' @param par The vector of parameters.
loglik <- function(p, y, X, par){
  # Parameters
  K <- length(p) # The number of mixture components
  B <- ncol(X) # The number of covariates
  beta <- matrix(par[1:(B * K)], nrow = B, ncol = K) # B*K betas
  phi <- par[((B * K + 1) : (B * K + K))] # K auto-regressive parameters
  kappa <- exp(par[(B * K + K + 1) : (B * K + K + K)]) # K concentrations
  
  eta <- 2 * atan(X %*% beta) %>% unname() # linear predictor
  r <- get_lagged_res(y, eta, burst, id) # lagged residuals
  kappa_t <- sqrt(kappa^2 + (phi * sin(r))^2) # time specific kappa
  mu_t <- eta + atan( (phi * sin (r)) / kappa_t) # time specific mu
  
  L <- sapply(
    seq_len(K),
    get_L,
    p = p, 
    y = y, 
    mu_t = mu_t, 
    kappa_t = kappa_t
    )
  
  if(nrow(r[!complete.cases(r),]) != nrow(L[!complete.cases(L),])){
    warning(
      "# of groups:", nrow(r[!complete.cases(r),]), ", ",
      "# of NAs:", nrow(L[!complete.cases(L),]))
  }
  
  L <- na_omit(L)
  
  list(L, 
       sum(log(rowSums(L))))
}

#' EM 
em <- function(y, X, K) {
  # Initial values
  par <- c(rnorm(K * ncol(X)), rnorm(K) ,rexp(K))
  p <- rexp(K)
  p <- p/sum(p)
  
  # EM-algorith
  # Initial likelihood 
  l <- vector()
  l[1] <- -Inf
  l[2] <- loglik(p, y, X, par)[[2]]
  
  it <- 2
  
  while (abs(l[it] - l[it - 1]) >= 1e-4){
    # E step
    f <- loglik(p, y, X, par)[[1]]
    w <- f/rowSums(f)
    
    # M-step 
    # pi_k explicit update
    p <- sapply(seq_along(p), function(k) mean(w[, k]))
    
    # Parameters computational update
    opt <- optim(
      par = par, 
      fn = m, 
      w = w, 
      y = y, 
      X = X,
      burst = burst, 
      id = id,
      method = "BFGS",
      control = list(maxit = 10000)
    )
    
    if (opt$convergence != 0) {# Convergence check
      warning("Iteration: ", it, "; Convergence: ", opt$convergence)
    }
    
    # Update iteration value
    par <- opt$par
    it <- it +1
    l[it] <- loglik(p, y, X, par)[[2]]
  }
  
  list(p=p, 
       par = par, 
       l = l)
}

#' ----------------------------------------------------------------------------
#' Data load and results

dir_data <- paste0(getwd(), "/data")
load(paste0(dir_data, "/elephants.RData"))

# Design matrix 
formula <- ~ distriv_std + elev_std + ndvi_std + seas + sex

data <- data %>%
  dplyr::select(all.vars(formula), "ta_", "burst_", "id") %>%
  tidyr::drop_na()

X <-  model.matrix(formula, data)
y <- data$ta_
burst <- data$burst_
id <- as.factor(data$id)


#' TODO
#' -[ ] Componente autoregressiva 
#'  -[ ] `get_lagged_res` vettoriale (più efficiente)
#'  -[ ] `get_q` e `get_loglik` esterni
#' -[ ] Stima per 3, 4 componenti 