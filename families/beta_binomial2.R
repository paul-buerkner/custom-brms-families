# This is the beta-binomial distribution in mean-precision parameterization
# Details are provided in
# https://paul-buerkner.github.io/brms/articles/brms_customfamilies.html

library(brms)

# helper functions for post-processing of the family
log_lik_beta_binomial2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}

posterior_epred_beta_binomial2 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  trials <- prep$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

# definition of the custom family
beta_binomial2 <- function(link = "logit", link_phi = "log") {
  custom_family(
    "beta_binomial2",
    dpars = c("mu", "phi"),
    links = c(link, link_phi),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "int",
    vars = "trials[n]",
    log_lik = log_lik_beta_binomial2,
    posterior_predict = posterior_predict_beta_binomial2,
    posterior_epred = posterior_epred__beta_binomial2
  )
}

# additionally required Stan code
stan_beta_binomial2 <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
