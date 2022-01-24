library(brms)

logit <- function(x) {
  log(x) - log1p(-x)
}

logistic <- function(x) {
  1 / (1 + exp(-x))
}

dlogitnormal <- function(x, md, sigma, log = FALSE){
  if (isTRUE(any(x <= 0 || x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(md <= 0 || md >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  if (log) {
    (- (log(sigma) + 0.5 * (log(2) + log(pi)))) +
    (- (log(x) + log1p(-x))) +
    (- (logit(x) - logit(md))^2) / (2 * (sigma^2))
  }
  else {
    exp(
      (- (log(sigma) + 0.5 * (log(2) + log(pi)))) +
      (- (log(x) + log1p(-x))) +
      (- (logit(x) - logit(md))^2) / (2 * (sigma^2))
    )
  }
}

rlogitnormal <- function(n, md, sigma) {
  if (isTRUE(any(md <= 0 || md >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    logistic(rnorm(n, logit(md), sigma))
  )
}

log_lik_logitnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  dlogitnormal(y, mu, sigma, log = TRUE)
}

posterior_predict_logitnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  rlogitnormal(1, mu, sigma)
}

posterior_epred_logitnormal <- function(prep) {
  stop("Due to the mean not having an analytical solution for the logit-normal
        distribution, posterior_epred is currently not supported.")
}

logitnormal <- custom_family(
  "logitnormal",
  dpars = c("mu", "sigma"),
  links = c("logit", "log"),
  lb = c(0, 0),
  ub = c(1, NA),
  type = "real",
  log_lik = log_lik_logitnormal,
  posterior_predict = posterior_predict_logitnormal,
  posterior_epred = posterior_epred_logitnormal
)

stan_logitnormal <- "
  real logitnormal_lpdf(real y, real mu, real sigma) {
    return log(1/(sigma * sqrt(2 * pi()))) + log(1/(y * (1-y))) +
           ((-(logit(y) - logit(mu))^2)/(2*(sigma^2)));
  }

  real logitnormal_rng(real mu, real sigma) {
    return inv_logit(normal_rng(logit(mu), sigma));
  }
"