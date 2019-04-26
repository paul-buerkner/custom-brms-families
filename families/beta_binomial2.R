# This is the beta-binomial distribution in mean-precision parameterization
# Details are provided in
# https://paul-buerkner.github.io/brms/articles/brms_customfamilies.html

library(brms)

# helper functions for post-processing of the family
log_lik_beta_binomial2 <- function(i, draws) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  N <- draws$data$trials[i]
  y <- draws$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, N)
}

predict_beta_binomial2 <- function(i, draws, ...) {
  mu <- draws$dpars$mu[, i]
  phi <- draws$dpars$phi
  N <- draws$data$trials[i]
  beta_binomial2_rng(mu, phi, N)
}

fitted_beta_binomial2 <- function(draws) {
  mu <- draws$dpars$mu
  trials <- draws$data$trials
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

# definition of the custom family
beta_binomial2 <- custom_family(
  "beta_binomial2", 
  dpars = c("mu", "phi"),
  links = c("logit", "log"), 
  lb = c(0, 0), 
  ub = c(1, NA),
  type = "int", vars = "trials[n]",
  log_lik = log_lik_beta_binomial2,
  predict = predict_beta_binomial2,
  fitted = fitted_beta_binomial2
)

# additionally required Stan code
stan_beta_binomial2 <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
