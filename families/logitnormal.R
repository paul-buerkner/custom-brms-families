library(brms)

logit <- function(x) {
  log(x) - log1p(-x)
}

logistic <- function(x) {
  1 / (1 + exp(-x))
}

dlogitnormal <- function(x, mu, sigma, log = FALSE){
  if (isTRUE(any(x <= 0 || x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  res <- (- (log(sigma) + 0.5 * (log(2) + log(pi)))) +
         (- (log(x) + log1p(-x))) +
         (- (logit(x) - mu)^2) / (2 * (sigma^2))
  if (log) {
    return(res)
  }
  else {
    return(exp(res))
  }
}

rlogitnormal <- function(n, mu, sigma) {
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    logistic(rnorm(n, mu, sigma))
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
  rlogitnormal(prep$ndraws, mu, sigma)
}

posterior_epred_logitnormal <- function(prep) {
  # https://doi.org/10.1080/03610926.2020.1752723 might solve this
  stop("Due to the mean not having an analytical solution for the logit-normal
        distribution, posterior_epred is currently not supported.")
}

logitnormal <- function(link = "identity", link_sigma = "log"){
  custom_family(
    "logitnormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_logitnormal,
    posterior_predict = posterior_predict_logitnormal,
    posterior_epred = posterior_epred_logitnormal
  )
}

stan_logitnormal <- "
  real logitnormal_lpdf(real y, real mu, real sigma) {
    return log(1/(sigma * sqrt(2 * pi()))) + log(1/(y * (1-y))) +
           ((-(logit(y) - logit(mu))^2)/(2*(sigma^2)));
  }

  real logitnormal_rng(real mu, real sigma) {
    return inv_logit(normal_rng(logit(mu), sigma));
  }
"
