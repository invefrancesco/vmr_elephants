setwd(rprojroot::find_rstudio_root_file())
source("code/_master_code.R")
#'
#' ### Introduzione 
#' 
#' La ragione per cui i dati fittati sono concentrati intorno alla media sembra essere che il modello ha una bassa capacità esplicativa. Essendo il valore del parametro stimato (in particolare quello della distanza dalla riva) basso i dati fittati non si allontanano dalla media. 
#' Faccio due ulteriori test qui per verificare che il problema sia quello
#' 
#' 1. Replico il modello con dati simulati. In teoria il grafico dei dati fittati dovrebbe avere lo stesso supporto dei dati reali.
#' 2. Replico il modello sui dati di un elefante con osservazioni più frequenti. Mi aspetto che se le osservazioni sono più frequenti le covariate siano più "capaci" di spiegare il movimento.
#' 
#' ### Primo test 
#| label: primo-test

# DATI RANDOM ----
set.seed(123)
dt <- tibble(
  x = rnorm(1000),
  y = sapply(x, function(i) circular::rvonmises(n = 1, mu = 2*atan(1.5 + 3.5 * i), kappa = 10)))

# OTTIMIZZAZIONE DELLA FUNZIONE 
results <- optim(
  par = c(0,0,log(2)),
  fn = log.lik.VM,
  data = dt,
  formula = y ~ x, 
  hessian = T
)

tibble(
  parametri = c("beta 0", "beta 1", "kappa"), 
  estimate = c(results$par[1:2], exp(results$par[3])) ,
  se = sqrt(diag(solve(results$hessian))),
  lower = estimate - qnorm(0.975) * se,
  upper = estimate + qnorm(0.975) * se
) %>% 
  kable(
    col.names = c("Parameter", "Estimate", "Std. Error", "95% CI Lower", "95% CI Upper"),
    align = "lcccc",
    format = "markdown", 
    digits = 3
  )
)

# DATI FITTATI 

dm <- mold(y~x, dt, 
           blueprint = default_formula_blueprint(intercept = TRUE))

# PARAMETERS
p <- ncol(design_matrix$predictors)
beta <- results$par[1:p]
kappa <- results$par[p+1]

# FUNCTION
eta <- as.vector(as.matrix(dm$predictors) %*% beta)

l <- tibble(
  y = as.vector(dm$outcomes$y),
  eta = eta
) %>% 
  mutate(
    mu = 2 * atan(eta),
    mu_std = ifelse(mu < 0, mu + 2*pi, mu)
  )

# Istogramma Fitted vs Observed 
ggplot(l, aes(x = y)) +
  geom_histogram(aes(y = ..density.., fill = "Observed"),
                 bins = 30, 
                 position = "identity", 
                 color = NA) +
  geom_histogram(aes(x = mu_std, y = ..density.., fill = "Fitted"), 
                 alpha = 0, 
                 bins = 30, 
                 position = "identity", 
                 color = "black") +
  scale_fill_manual(
    values = c("Observed" = "red", "Fitted" = "black"),
    name = NULL
  ) +
  labs(
    title = "Fitted vs Observed",
    x = "Angle (radians)",
    y = "Density"
  ) +
  theme_test()

