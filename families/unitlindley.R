# suggested in doi:10.1080/02664763.2018.1511774 for data on the unit interval

library(brms)
library(lamW) # for the quantile function in posterior_predict

# helper functions for post-processing of the family
log_lik_unitlindley <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  y <- prep$data$Y[i]
  2*log(1-mu)-log(mu)-3*log(1-y)-y*(1-mu)/(mu*(1-y))
}

posterior_epred_unitlindley <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

posterior_predict_unitlindley <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  lambert <- lambertWm1((runif(1, 0, 1)-1)/mu*exp(-1/mu))
  return((1/mu+lambert)/(1+lambert))
}

# definition of the custom family

unitlindley <- custom_family("unitlindley", 
                             dpars = "mu", 
                             links = "logit", 
                             lb = 0, 
                             ub = 1, 
                             type = "real", 
                             log_lik = log_lik_unitlindley, 
                             posterior_epred = posterior_epred_unitlindley,
                             posterior_predict = posterior_predict_unitlindley)
# additionally required Stan code

stan_funs <- "
real unitlindley_lpdf(real y, real mu){
  return 2*log(1-mu)-log(mu)-3*log(1-y)-y*(1-mu)/(mu*(1-y));
}
"

# template
# stanvars = stanvar(scode = stan_funs, block = "functions")
# brm(..., stanvars = stanvars, family = unitlindley)
