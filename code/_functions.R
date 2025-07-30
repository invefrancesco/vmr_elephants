# regressione Von Mises ----
log.lik.VM <- function(par, data, formula, response){
  # Dati
  X <- model.matrix(formula, data = data)
  y <- data[[response]]
  
  # Parametri
  p <- ncol(X)
  beta <- par[1:p]
  kappa <- exp(par[p+1])
  
  # Elimina righe con NA
  valid <- complete.cases(y, X)
  y <- y[valid]
  X <- X[valid, , drop = FALSE]
  
  # Link function
  eta <- X %*% beta
  mu <- 2 * atan(eta)
  
  # Log-likelihood function
  l <- kappa * cos(y - mu) - log(besselI(kappa, nu = 0))
  return(-sum(l))
}

# Circular Autocorrelation Function ----
acf_circular <- function(lag, x){
  res_lag <- x %>%  
    mutate(t2_lag = t2_ + hours(lag * 4))
  
  res_lag <- res_lag %>% 
    inner_join(
      residuals,
      by = c("t2_lag" = "t2_")
    )
  
  n_match <- nrow(res_lag)
  
  cor <-  circular::cor.circular(res_lag$res.x, res_lag$res.y, test = T)
  tibble(
    lag = lag,
    n_match = n_match,
    acf = cor[[1]],
    statistic = cor[[2]],
    p.value = cor[[3]]
  )
}

# log-likelihood AR(1) ----
log.lik.VM.ar <- function(par, data, formula, response, burst) {
  # dati ----
  X <- model.matrix(formula, data)
  y <- data[[response]]
  burst <- data[[burst]]
  
  # parametri ----
  beta.0 <- atan2(sin(par[1]), cos(par[1]))
  beta <- par[2:ncol(X)]
  phi <- par[ncol(X) + 1]
  kappa <- exp(par[ncol(X) + 2])
  
  # log-likelihood ----
  eta <- X[,-1] %*% beta
  l_data <- tibble(
    burst = burst,
    y = y,
    eta = as.vector(eta)
  ) %>% 
    group_by(burst) %>% 
    mutate(
      res = y - 2 * atan(eta),
      kappa_t = sqrt(kappa^2 + (phi * sin(lag(res)))^2),
      mu_t =  beta.0 + 2 * atan(eta) + atan(phi * sin(lag(res)) / kappa_t),
      l = kappa_t * cos(y - mu_t) - log(besselI(kappa_t, nu = 0))
    ) %>% 
    slice(3:n())
  
  return(-sum(l_data$l))
}