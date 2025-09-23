#| setup
source("code/_master_code.R")
#' 
#' # GON34
#' 
#' ## EDA
#' 
#| EDA
# data 
GON34 <- read_csv(paste0(dir_data, "/GON34.csv"))
# turn angle
ggplot(data = GON34, aes(x = ta_, y = ..density..)) +
  geom_histogram(color = "black", fill = NA) +
  geom_density() +
  facet_wrap(~ seas) +
  theme_test() +
  xlab("Turn angle") +
  ylab("Density")
#| fit 
fit <- optim(
  par = c(rnorm(4), log(2)),
  fn = log.lik.VM, 
  data = GON34, 
  formula = ta_ ~ distriv_std + elev_std + ndvi_std,
  control = list(maxit = 10000),
  hessian = T)
# risultati ----
#     tabella ----
tibble(
  parametri = c("Intercept", "Distance from water", "Elevation", "NDVI Index", "Kappa"),
  estimate = c(fit$par[1:4], exp(fit$par[5])), 
  se = sqrt(diag(solve(fit$hessian))),
  lower = estimate - qnorm(0.975) * se,
  upper = estimate + qnorm(0.975) * se,
  W = estimate / se, 
  p_value = 2 * (1 - pnorm(abs(W)))
) %>% 
  kable(
    col.names = c("Parameter", "Estimate", "Std. Error", "95% CI Lower", "95% CI Upper", "Wald test", "p-value"),
    align = "lcccc",
    format = "markdown", 
    
  )
#     grafico dei coefficienti ----
tibble(
  parameter = factor(c("Intercept", "Distance from water", "Elevation", "NDVI Index", "Kappa"),
                     levels = c("Intercept", "Distance from water", "Elevation", "NDVI Index", "Kappa")),
  estimate = c(fit$par[1:4], exp(fit$par[5])), 
  se = sqrt(diag(solve(fit$hessian))),
  lower = estimate - qnorm(0.975) * se,
  upper = estimate + qnorm(0.975) * se
) %>%
  ggplot(aes(estimate, parameter)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper)) +
  geom_vline(xintercept = 0, lty = 2, color = "red") +
  labs(
    x = "Estimate with conf. intervals"
  ) +
  theme_test()