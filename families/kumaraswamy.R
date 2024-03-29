library(brms)

dkumaraswamy <- function(x, mu, p, log = FALSE) {
  if (isTRUE(any(x <= 0 || x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(mu <= 0 || mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  res <- log(p) +
         log(log(2)) -
         log(- (log1p(-mu^p))) +
         (p - 1) * log(x) +
         ((- (log(2) / log1p(-mu^p))) - 1) * log1p(-x^p)

  if (log) {
    return(res)
  } else {
    return(exp(res))
  }
}

rkumaraswamy <- function(n, mu, p) {
  if (isTRUE(any(mu <= 0 || mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  q <- - (log(2) / log1p(-mu^p))
  return(
    (1 - (1 - runif(n, min = 0, max = 1))^(1 / q))^(1 / p)
  )
}

log_lik_kumaraswamy <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  y <- prep$data$Y[i]
  dkumaraswamy(y, mu, p, log = TRUE)
}

posterior_predict_kumaraswamy <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  rkumaraswamy(prep$ndraws, mu, p)
}

posterior_epred_kumaraswamy <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  p <- brms::get_dpar(prep, "p")
  q <- - (log(2) / log1p(-mu^p))
  q * beta((1 + 1 / p), q)
}

kumaraswamy <- function(link = "logit", link_p = "log"){
  custom_family(
    "kumaraswamy",
    dpars = c("mu", "p"),
    links = c(link, link_p),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_kumaraswamy,
    posterior_predict = posterior_predict_kumaraswamy,
    posterior_epred = posterior_epred_kumaraswamy
  )
}

stan_kumaraswamy <- "
  real kumaraswamy_lpdf(real y, real mu, real p) {
    return  (log(p) + log(log(2)) - log(-(log1m(mu^p))) + (p-1) * log(y) +
            ((-(log(2)/log1m(mu^p)))-1) * log1m(y^p));
  }

  real kumaraswamy_rng(real mu, real p) {
    return ((1-(1-uniform_rng(0, 1))^(1/(-(log(2)/log1m(mu^p)))))^(1/p));
  }
"
