#' # Script for MLE with fixed effect -----------------------------------------
#' 
#' Load packages, functions and data
source("code/_master_code.R")
load(paste0(dir_data, "/elephants.RData"))

#' Optimization of the baseline function, assuming iid angles -----------------
formula <- ~ sl_ + distriv_std + ndvi_std + elev_std + sex + seas

fit <- optim(
  fn = log.lik.VM.ar,
  par = c(rnorm(9), log(2)),
  data = data,
  formula = formula, 
  response = "ta_",
  burst = "burst_", 
  method = "L-BFGS-B",
  control = list(maxit = 1000),
  hessian = TRUE,
)

#' Visualization of coefficients ----------------------------------------------
coeff <- tibble(
  parameter = c(
    "Intercept",
    "Step Length",
    "Distance from water",
    "NDVI Index",
    "Elevation",
    "Sex = Male",
    "Season = HD",
    "Season = HW",
    "Phi",
    "Kappa"
  ),
  estimate = c(fit$par[1:9], exp(fit$par[10])),
  se = sqrt(diag(solve(fit$hessian))),
  lower = estimate - qnorm(0.975) * se,
  upper = estimate + qnorm(0.975) * se,
  W = estimate / se,
  p_value = 2 * (1 - pnorm(abs(W)))
)

kable(
  coeff,
  col.names = c(
    "Parameter",
    "Estimate",
    "Std. Error",
    "95% CI Lower",
    "95% CI Upper",
    "Wald test",
    "p-value"
  ),
  align = "lcccc",
  format = "markdown"
)

#' Optimization of the function with Fixed Effect -----------------------------
#' Complete cases
formula <- ta_ ~ sl_ + distriv_std + ndvi_std + elev_std + sex + seas + id
data_model <- data %>%# nolint
  dplyr::select(all.vars(formula), "ta_", "burst_", "id") %>%# nolint
  tidyr::drop_na()

#' Design matrix
mm <- model.matrix(formula, data_model)

#' Formula
formula2 <- ta_ ~ sl_ + distriv_std + ndvi_std + elev_std + sex + seas + id

init = Sys.time()
fit2 <- optim(
  fn = log.lik.VM.ar.fe,
  par = c(rnorm(9), rnorm(23), log(2)),
  data = data_model,
  mm = mm, 
  response = "ta_",
  burst = "burst_", 
  id = "id", 
  method = "L-BFGS-B",
  hessian = FALSE, 
  control = list(maxit = 10, trace = 1)
)
Sys.time() - init
log.lik.VM.ar.fe(par = c(rnorm(9), rnorm(23), log(2)),
                 data = data_model, 
                 mm = mm,
                 response = "ta_", 
                 burst = "burst_",
                 id = "id",) 

# NOTE: The optim function takes >1h to converge