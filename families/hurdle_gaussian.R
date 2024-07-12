# https://discourse.mc-stan.org/t/custom-gaussian-hurdle-family-not-quite-working-in-brms/21028
library(brms)

# Create a custom family that is logit if y = 0, normal/gaussian if not
hurdle_gaussian <-
    custom_family("hurdle_gaussian",
        dpars = c("mu", "sigma", "hu"),
        links = c("identity", "log", "logit"),
        lb = c(NA, NA, 0),
        type = "real"
    )

# helper functions for post processing
dhurdle_gaussian <- function(x, mu = 0, sigma = 1, hu = 0, log = FALSE) {
    # based on https://www.m-flynn.com/posts/hurdle-lognormal-densities-ii/lognormal-densities-ii
    if (x > 0) {
        value <- dnorm(x = x, mean = mu, sd = sigma, log = log)
        return(value)
    } else {
        value <- hu
        return(value)
    }
}

stan_funs <- "
  real hurdle_gaussian_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             normal_lpdf(y | mu, sigma);
    }
  }
"

stanvars <- stanvar(scode = stan_funs, block = "functions")

posterior_predict_hurdle_gaussian <- function(i, prep, ...) {
    theta <- prep$dpars$hu[, i]
    mu <- prep$dpars$mu[, i]
    sigma <- prep$dpars$sigma
    ndraws <- prep$ndraws

    hu <- runif(n = ndraws, min = 0, max = 1)
    ifelse(hu < theta, 0, rnorm(ndraws, mu, sigma))
}

posterior_epred_hurdle_gaussian <- function(prep) {
    with(prep$dpars, mu * (1 - hu))
}


log_lik_weight <- function(x, i, prep) {
    weight <- prep$data$weights[i]
    if (!is.null(weight)) {
        x <- x * weight
    }
    x
}

log_lik_hurdle_gaussian <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    hu <- brms::get_dpar(prep, "hu", i = i)
    y <- prep$data$Y[i]
    dhurdle_gaussian(x = y, mu = mu, sigma = sigma, hu = hu, log = TRUE)
}

if (FALSE) {
    # Example
    library(data.table)
    library(ggplot2)

    example_data <- data.table("outcome" = rnorm(1000, 15, 3))

    example_data[, x := rnorm(1000, 5, 1)]
    example_data[, outcome := outcome * rbinom(1000, 1, 0.9) * (0.25 * x)]

    # Normal around 20ish, but lots of 0s
    ggplot(example_data, aes(x = outcome)) +
        geom_histogram(binwidth = 1, boundary = 0, color = "white")

    prop_zero <- function(x) mean(x == 0)
    zeroes_in_data <- prop_zero(example_data$outcome) # ~ 11% zeros
    logit_scaled(zeroes_in_data) # On logit scale, around -2

    model_fit <- brm(
        bf(
            outcome ~ x,
            hu ~ 1
        ),
        data = example_data,
        family = hurdle_gaussian,
        stanvars = stanvars,
        chains = 4,
        iter = 1000,
        warmup = 500,
        cores = 4,
        backend = "cmdstanr"
    )

    summary(model_fit)

    pp_check(model_fit, ndraws = 500)

    conditional_effects(model_fit)
}
