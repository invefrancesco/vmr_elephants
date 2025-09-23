# Regressione Von Mises ----
log.lik.VM <- function(par, data, formula) {
  # Preparazione dati 
  mm <- model.matrix(formula, data)
  mf <- model.frame(formula, data)
  
  # Parametri
  p     <- ncol(mm)      # numero di regressori
  beta  <- par[1:p]                 # coefficienti
  kappa <- exp(par[p + 1])          # forziamo positività con exp()
  
  # Predittore lineare
  eta <- mm %*% beta
  
  # log-likelihood
  y   <- mf[1]
  mu  <- 2 * atan(eta)
  
  # log-likelihood (usiamo expon.scaled per stabilità numerica)
  l <- kappa * cos(y - mu) - besselI(kappa, nu = 0)
  
  # Ritorniamo il negativo (per ottimizzazione)
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
  beta <- par[1:ncol(X)]
  phi <- par[ncol(X) + 1]
  kappa <- exp(par[ncol(X) + 2])
  
  # log-likelihood ----
  eta <- X %*% beta
  l_data <- tibble(
    burst = burst,
    y = y,
    eta = as.vector(eta)
  ) %>% 
    group_by(burst) %>% 
    mutate(
      res = y - 2 * atan(eta),
      kappa_t = sqrt(kappa^2 + (phi * sin(lag(res)))^2),
      mu_t =  2 * atan(eta) + atan(phi * sin(lag(res)) / kappa_t),
      l = kappa_t * cos(y - mu_t) - log(besselI(kappa_t, nu = 0))
    ) %>% 
    slice(3:n())
  
  return(-sum(l_data$l))
}

# log-likelihood ar alternativa 
log.lik.VM.ar.alt <- function(par, data, formula, response, burst, L) {
  # dati ----
  X <- model.matrix(formula, data)
  y <- data[[response]]
  burst <- data[[burst]]
  y_lag <- as_tibble(set_names(
    map(1:L, ~ lag(y, .x)),
    paste0("lag_", 1:L))) %>% 
    as.matrix()
  
  # parametri ----
  beta <- par[1:ncol(X)]
  phi <- par[(ncol(X) +1):(ncol(X) + L)]
  kappa <- exp(par[ncol(X) + L + 1])
  
  # log-likelihood ----
  eta <- X %*% beta
  ar_term <- y_lag %*% phi
  l_data <- tibble(
    burst = burst,
    y = y,
    eta = as.vector(eta),
    ar_term = as.vector(ar_term)
  ) %>% 
    group_by(burst) %>% 
    mutate(
      mu_t =  2 * atan(eta) + ar_term, 
      l = kappa * cos(y - mu_t) - log(besselI(kappa, nu = 0))
    ) %>% 
    slice((L+2):n())
  
  return(-sum(l_data$l))
}

# log-likelihood con FE
llVM.ar.fe <- function(par, data, formula, burst, id) {
  # data ----
  dm <- mold(formula, data,
             blueprint = default_formula_blueprint(intercept = T))
  
  # parameters ----
  beta <- par[1:ncol(dm$predictors)]
  phi <- par[ncol(dm$predictors) + 1] 
  kappa <- exp(par[ncol(dm$predictors) + 2])
  
  # log-likelihood ---
  eta <- as.matrix(dm$predictors) %*% beta # linear predictor
  ar_term <- as.matrix(lag(dm$outcomes)) %*% phi # autoregressive term 
  l_data <- tibble(
    id = data[[id]],
    burst = data[[burst]], 
    y = dm$outcomes[,1], 
    eta = as.vector(eta), 
    ar_term = as.vector(ar_term)
  ) %>% 
    group_by(id, burst) %>% # the likelihood is estimated inside id and burst 
    mutate(
      mu_t =  2 * atan(eta) + ar_term, 
      l = kappa * cos(y - mu_t) - log(besselI(kappa, nu = 0))
    ) %>% 
    slice(3:n()) # first two terms are (for id and burst) NA
  
  # response 
  return(-sum(l_data$l))
}