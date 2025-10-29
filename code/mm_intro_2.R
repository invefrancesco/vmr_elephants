# This is the second part of mm_intro 
# The previous analysis to a more complex setting:
#   1. Our data are now generated from a mixture of 3 components instead of two. 
#   2. The number of components assumed is unknown: AIC and BIC are used in
#   order to identify the best model.
#   3. The components are Von Mises distributions 
source("code/_master_code.R")
set.seed(12345)

# Generate a sample from three Von Mises --------------------------------------
u <- sample(1:3, 5000, replace = TRUE, prob = c(0.40, 0.35, 0.25))# categories
x <- vector()# empty vector

for(i in seq_along(u)) {
  if (u[i] == 1){
    x[i] <- circular::rvonmises(1, mu = pi, kappa = 10)
  } else if (u[i] == 2) {
    x[i] <- circular::rvonmises(1, mu = 0, kappa = 10) 
  } else {
    x[i] <-  circular::rvonmises(1, mu = pi/2, kappa = 10)
  }
}

# Visualization ---------------------------------------------------------------
x <- circular(x, 
              type = "angles", 
              units = "radians",
              modulo = "2pi",
              rotation = "counter")

ggplot(tibble(angle = x * 180 / pi), aes(x = angle, y = after_stat(density))) +
  geom_histogram(bins = 40, color = "black", fill = "gray") + 
  geom_density(color = "red", lty = 2) +
  coord_polar(start = 3*pi/2, direction = -1) +
  scale_x_continuous(
    limits = c(0, 360),
    breaks = c(0, 90, 180, 270),
    labels = c("0", "π/2", "π", "3π/2")
  ) +
  theme_bw() +
  labs(x = NULL, y = NULL)

# Initialize parameters -------------------------------------------------------
# weights
pi <- rexp(2)
pi <- pi/sum(pi)
# mu
mu <- rnorm(2)
kappa <- rep(10, 2) #start with k given

# E-step ----------------------------------------------------------------------
# Posterior probabilities w[i, k] are given by joint/marginal distribution
# evaluated in x[i] given the initial parameters pi[k], mu[k], kappa[k]
f <- sapply(seq_along(pi), function(k) pi[k] * dvonmises(x, mu[k], kappa[k]))
w <- f/rowSums(f)

# w is a matrix n * k. As a check, we can verify that sum(w) = n, since each 
# row of the weights mus sum to 1
sum(w)

# M-step -----------------------------------------------------------------------
# The update for pi[k] is always given by 1/n * sum(w[ ,k]) = mean(w[ ,k])
pi <- sapply(seq_along(pi), function(k) mean(w[, k]))

# The update for mu[k] is given by:
# sum(w[i, k] * sin(x[i])) / sum(w[i, k] * cos(x[i]))
mu <- sapply(seq_along(mu), 
             function(k) sum(w[ ,k] * sin(x))/ sum(w[ ,k] * cos(x))
)

# The update for the parameters could also be found by the numerical 
# optimization of the second part of the complete likelihood (given w).
# In this case it is the sum over all the i and k of the function: 
# w[i, k] * dvonmises(x[i], mu[k], kappa[k], log = TRUE)

m <- function(w, x, par){
  mu <- par[1:ncol(w)]
  kappa <- exp(par[(ncol(w)+1):(ncol(w)*2)])
  
  Q <- sapply(seq_along(mu), 
              function(k) w[ ,k] * dvonmises(x, mu[k], kappa[k], log = TRUE)
  )
  
  -sum(Q)
}

# Log-likelihod ----------------------------------------------------------------
# In order to verify the convergence of the algorithm we store the value of 
# the log-likelihood function at each iteration

loglik <- function(x, mu, kappa, pi){
  # Compute the value of each component of the log-likelihood
  # It must be a matrix n * k
  l <- sapply(seq_along(pi), function(k) pi[k] * dvonmises(x, mu[k], kappa[k]))
  
  # Compute the sum of the components
  sum(log(rowSums(l)))
}

l <- vector()
l[1] <- -Inf
l[2] <- loglik(x, mu, kappa, pi)

# EM algorithm ----------------------------------------------------------------
for(t in 3:30) {
  # E step
  f <- sapply(seq_along(pi), function(k) pi[k] * dvonmises(x, mu[k], kappa[k]))
  w <- f/rowSums(f)
  
  # M-step
  pi <- sapply(seq_along(pi), function(k) mean(w[, k]))
  mu <- sapply(seq_along(mu),
               function(k) atan2(sum(w[ ,k] * sin(x)), sum(w[ ,k] * cos(x)))
  )
  
  l[t] <- loglik(x, mu, kappa, pi)
}

plot(l)

# EM algorithm with optim -----------------------------------------------------
rm(list = c("f", "w", "i", "kappa", "l", "mu", "pi", "t", "u"))

# Initial parameters (2 components)
pi <- rexp(2)
pi <- pi/sum(pi) 
mu <- rnorm(2)
kappa <- abs(rnorm(2)) #start with k given

l <- vector()
l[1] <- -Inf
l[2] <- loglik(x, mu, kappa, pi)

it <- 2
while(abs(l[it] - l[it - 1]) >= 1e-6) {
  # E step
  f <- sapply(seq_along(pi), function(k) pi[k] * dvonmises(x, mu[k], kappa[k]))
  w <- f/rowSums(f)
  
  # M-step 
  pi <- sapply(seq_along(pi), function(k) mean(w[, k]))
  
  opt <- optim(
    fn = m,
    par = c(mu, log(kappa)),
    w = w,
    x = x
  )
  
  mu <- opt$par[1:length(pi)]
  kappa <- exp(opt$par[(length(pi) + 1) : (length(pi) * 2)])
  
  it <- it +1
  l[it] <- loglik(x, mu, kappa, pi)
}

plot(l)

em_results <- list()
em_results[["2"]] <- list(
  K = 2,
  pi = pi,
  mu = mu,
  kappa = kappa,
  loglik = l[it],
  loglik_trace = l
)

# 3 components 
rm(list = c("f", "w", "i", "kappa", "l", "mu", "pi", "t", "u"))

# Initial parameters (2 components)
pi <- rexp(3)
pi <- pi/sum(pi) 
mu <- rnorm(3)
kappa <- abs(rnorm(3)) #start with k given

l <- vector()
l[1] <- -Inf
l[2] <- loglik(x, mu, kappa, pi)

it <- 2
while(abs(l[it] - l[it - 1]) >= 1e-6) {
  # E step
  f <- sapply(seq_along(pi), function(k) pi[k] * dvonmises(x, mu[k], kappa[k]))
  w <- f/rowSums(f)
  
  # M-step 
  pi <- sapply(seq_along(pi), function(k) mean(w[, k]))
  
  opt <- optim(
    fn = m,
    par = c(mu, log(kappa)),
    w = w,
    x = x
  )
  
  mu <- opt$par[1:length(pi)]
  kappa <- exp(opt$par[(length(pi) + 1) : (length(pi) * 2)])
  
  it <- it +1
  l[it] <- loglik(x, mu, kappa, pi)
}
plot(l)

em_results[["3"]] <- list(
  K = 3,
  pi = pi,
  mu = mu,
  kappa = kappa,
  loglik = l[it],
  loglik_trace = l
)

# 4 components


# TODO ------------------------------------------------------------------------
# - [ ] AIC e BIC
# - [ ] Bootstrap SE