# Von Mises log-likelihood ----------------------------------------------------
log.lik.VM <- function(par, data, formula, response) {
  # Complete cases
  data <- data %>% # nolint
    dplyr::select(all.vars(formula), dplyr::all_of(response)) %>%
    tidyr::drop_na()
  # Preparazione dati
  mm <- model.matrix(formula, data)
  y <- data[[response]]

  # Parametri
  beta <- par[seq_len(ncol(mm))] # coefficienti
  kappa <- exp(par[ncol(mm) + 1]) # forza la positività di kappa

  # Predittore lineare
  eta <- mm %*% beta

  # log-likelihood
  mu <- 2 * atan(eta)
  l <- kappa * cos(y - mu) - besselI(kappa, nu = 0)

  -sum(l)
}

# Von Mises log-likelihood with autoregressive component AR(1) ----------------
log.lik.VM.ar <- function(par, data, formula, response, burst) {
  # Complete cases
  data <- data %>% # nolint
    dplyr::select(all.vars(formula), dplyr::all_of(c(response, burst))) %>%# nolint
    tidyr::drop_na()

  # Design matrix
  mm <- model.matrix(formula, data)
  y <- data[[response]]
  burst <- data[[burst]]

  # Parameters
  beta <- par[seq_len(ncol(mm))]
  phi <- par[ncol(mm) + 1]
  kappa <- exp(par[ncol(mm) + 2]) # forza la positività di kappa

  # log-likelihood ----
  eta <- drop(mm %*% beta)
  l_data <- tibble(
    burst = burst,
    y = y,
    eta = eta
  ) %>%
    group_by(burst) %>%
    mutate(
      res = y - 2 * atan(eta),
      kappa_t = sqrt(kappa^2 + (phi * sin(lag(res)))^2),
      mu_t = 2 * atan(eta) + atan(phi * sin(lag(res)) / kappa_t),
      l = kappa_t * cos(y - mu_t) - log(besselI(kappa_t, nu = 0))
    ) %>%
    slice(2:n())

  -sum(l_data$l)
}

# Von Mises log-likelihood with autoregressive component AR(1) ----------------
# In this version, a id column is added in order to estimate FE
log.lik.VM.ar.fe <- function(par, data, formula, response, burst, id) {
  # Complete cases
  data <- data %>%# nolint
    dplyr::select(all.vars(formula), dplyr::all_of(c(response, burst))) %>%# nolint
    tidyr::drop_na()
  
  # Design matrix
  mm <- model.matrix(formula, data)
  y <- data[[response]]
  burst <- data[[burst]]
  id <- data[[id]]
  
  # Parameters
  beta <- par[seq_len(ncol(mm))]
  phi <- par[ncol(mm) + 1]
  kappa <- exp(par[ncol(mm) + 2]) # forza la positività di kappa
  
  # log-likelihood ----
  eta <- drop(mm %*% beta)
  l_data <- tibble(
    burst = burst,
    id = id,
    y = y,
    eta = eta
  ) %>%
    group_by(id, burst) %>%
    mutate(
      res = y - 2 * atan(eta),
      kappa_t = sqrt(kappa^2 + (phi * sin(lag(res)))^2),
      mu_t = 2 * atan(eta) + atan(phi * sin(lag(res)) / kappa_t),
      l = kappa_t * cos(y - mu_t) - log(besselI(kappa_t, nu = 0))
    ) %>%
    slice(2:n())
  
  -sum(l_data$l)
}
