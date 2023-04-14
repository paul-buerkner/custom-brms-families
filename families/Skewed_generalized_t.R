library(brms)

## Helper functions for post-processing

rskew_generalized_t <- function(n, mu, sigma, lambda, p, q) {
  # Generate random samples from a standard skew generalized t distribution
  u <- runif(n)
  v <- qbeta(u, 1/p, q)  
  z <- qnorm(u * pnorm((1 - lambda) * sqrt(v), lower.tail = FALSE))
  r <- sqrt(v) * (z + lambda * abs(z))
  
  # Scale and shift the random samples
  return(mu + sigma * r)
}


# Helper function for lpdf
dskew_generalized_t <- function(x, mu, sigma, lambda, p, q, log = FALSE) {
  stan_fit <- rstan::stan(model_code = stan_code)
  res <- rstan::expose_stan_functions(stan_fit)
  log_prob <- res$skew_generalized_t_lpdf(x, mu, sigma, lambda, p, q)
  
  if (log) {
    return(log_prob)
  } else {
    return(exp(log_prob))
  }
}

# Log likelihood function
log_lik_skew_generalized_t <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  lambda <- brms::get_dpar(prep, "lambda", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  q <- brms::get_dpar(prep, "q", i = i)
  y <- prep$data$Y[i]
  
  return(dskew_generalized_t(y, mu, sigma, lambda, p, q, log = TRUE))
}


# Posterior prediction function
posterior_predict_skew_generalized_t <- function(i, prep, ..., ntrys = 5) {
  # Add the internal functions get_se, add_sigma_se and rcontinuous
  # get_se, add_sigma_se and rcontinuous copied and adapted from brms repository 
  get_se <- function(prep, i = NULL) {
    stopifnot(is.brmsprep(prep))
    se <- as.vector(prep$data[["se"]])
    if (!is.null(se)) {
      if (!is.null(i)) {
        se <- se[i]
      }
      if (length(se) > 1L) {
        dim <- c(prep$ndraws, length(se))
        se <- data2draws(se, dim = dim)
      }
    } else {
      se <- 0
    }
    se
  }
  
  add_sigma_se <- function(sigma, prep, i = NULL) {
    if ("se" %in% names(prep$data)) {
      se <- get_se(prep, i = i)
      sigma <- sqrt(se^2 + sigma^2)
    }
    sigma
  }
  
  rcontinuous <- function(n, dist_name, dist, ..., lb = NULL, ub = NULL, ntrys = 5) {
    args <- list(...)
    if (is.null(lb) && is.null(ub)) {
      # sample as usual
      rdist <- paste0("r", dist_name)
      out <- do_call(rdist, c(list(n), args))
    } else {
      # sample from truncated distribution
      pdist <- paste0("p", dist_name)
      qdist <- paste0("q", dist_name)
      if (!exists(pdist, mode = "function") || !exists(qdist, mode = "function")) {
        # use rejection sampling as CDF or quantile function are not available
        out <- rdiscrete(n, dist, ..., lb = lb, ub = ub, ntrys = ntrys)
      } else {
        if (is.null(lb)) lb <- -Inf
        if (is.null(ub)) ub <- Inf
        plb <- do_call(pdist, c(list(lb), args))
        pub <- do_call(pdist, c(list(ub), args))
        out <- runif(n, min = plb, max = pub)
        out <- do_call(qdist, c(list(out), args))
        # infinite values may be caused by numerical imprecision
        out[out %in% c(-Inf, Inf)] <- NA
      }
    }
    out
  }
  
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  sigma <- add_sigma_se(sigma, prep, i = i) # handle sigma se
  lambda <- brms::get_dpar(prep, "lambda", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  q <- brms::get_dpar(prep, "q", i = i)
  
  rcontinuous(
    n = prep$ndraws,
    dist_name = "skew_generalized_t", # pass the distribution name separately
    dist = function(n, mu, sigma, lambda, p, q) rskew_generalized_t(n, mu, sigma, lambda, p, q),
    mu = mu,
    sigma = sigma,
    lambda = lambda,
    p = p,
    q = q,
    lb = prep$data$lb[i],
    ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

# Posterior expected prediction function
posterior_epred_skew_generalized_t <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  lambda <- brms::get_dpar(prep, "lambda")
  p <- brms::get_dpar(prep, "p")
  q <- brms::get_dpar(prep, "q")
  
  ndraws <- prep$ndraws
  ns <- nrow(mu)
  
  # Use apply() to generate the epred matrix
  epred <- t(apply(cbind(mu, sigma, lambda, p, q), 1, function(x) {
    rskew_generalized_t(ndraws, x[1], x[2], x[3], x[4], x[5])
  }))
  
  
  return(epred)
}

# Custom family definition
skew_generalized_t <- custom_family("skew_generalized_t",
                                    dpars = c("mu", "sigma", "lambda", "p", "q"),
                                    links = c("identity", "log", "identity", "identity", "identity"),
                                    type = "real",
                                    posterior_predict = posterior_predict_skew_generalized_t,
                                    posterior_epred = posterior_epred_skew_generalized_t)


# Custom family functions
skew_generalized_t$functions <- list(log_lik = "log_lik_skew_generalized_t")


##Initial values for skew_generalized_t parameters

init_fun <- function() {
  list(mu = rnorm(1, 0, 1),
       sigma = runif(1, 0.1, 1),
       lambda = runif(1, -0.99, 0.99),
       p = runif(1, 2, 10),
       q = runif(1, 2, 10))
}


# Required Stan code
# Stan code provided by Sean Pinkney
# https://github.com/spinkney/helpful_stan_functions/blob/main/functions/distribution/skew_generalized_t.stanfunctions

stan_skew_generalized_t <- '
real variance_adjusted_sgt(real sigma, real lambda, real p, real q) {
    if (p * q <= 2) 
      reject("p * q must be > 2 found p * q = ", p * q);
    
    if (is_inf(q)) 
      return sigma
             * inv_sqrt((pi() * (1 + 3 * lambda ^ 2) * tgamma(3.0 / p)
                         - 16 ^ (1.0 / p) * lambda ^ 2 * (tgamma(1.0 / 2 + 1.0 / p)) ^ 2 * tgamma(1.0 / p))
                        / (pi() * tgamma(1.0 / p)));
    
    return sigma
           / (q ^ (1.0 / p)
              * sqrt((3 * lambda ^ 2 + 1) * (beta(3.0 / p, q - 2.0 / p) / beta(1.0 / p, q))
                     - 4 * lambda ^ 2 * (beta(2.0 / p, q - 1.0 / p) / beta(1.0 / p, q)) ^ 2));
  }
  
  vector mean_centered_sgt(vector x, real sigma, real lambda, real p, real q) {
    if (p * q <= 1) 
      reject("p * q must be > 1 found p * q = ", p * q);
    
    if (is_inf(q)) 
      return x + (2 ^ (2.0 / p) * sigma * lambda * tgamma(1.0 / 2 + 1.0 / p)) / sqrt(pi());
    
    return x + (2 * sigma * lambda * q ^ (1.0 / p) * beta(2 / p, q - 1.0 / p)) / beta(1.0 / p, q);
  }
  
  real mean_centered_sgt(real x, real sigma, real lambda, real p, real q) {
    if (p * q <= 1) 
      reject("p * q must be > 1 found p * q = ", p * q);
    
    if (is_inf(q)) 
      return x + (2 ^ (2.0 / p) * sigma * lambda * tgamma(1.0 / 2 + 1.0 / p)) / sqrt(pi());
    
    return x + (2 * sigma * lambda * q ^ (1.0 / p) * beta(2 / p, q - 1.0 / p)) / beta(1.0 / p, q);
  }
  
  real mean_centered_sgt(real x, real sigma, real lambda, real q) {
    if (q <= 1) 
      reject("q must be > 1 found q = ", q);
    
    if (is_inf(q)) 
      return x + (4 * sigma * lambda * tgamma(1.0 / 2 + 1.0)) / sqrt(pi());
    
    return x + (2 * sigma * lambda * q * beta(2, q - 1.0)) / beta(1.0, q);
  }
  
  real skew_generalized_t_lpdf(real x, real mu, real sigma, real lambda, real p, real q) {
    if (sigma <= 0) 
      reject("sigma must be > 0 found sigma = ", sigma);
    
    if (lambda >= 1 || lambda <= -1) 
      reject("lambda must be between (-1, 1) found lambda = ", lambda);
    
    if (p <= 0) 
      reject("p must be > 0 found p = ", p);
    
    if (q <= 0) 
      reject("q must be > 0 found q = ", q);
    
    int N = 1;
    real out = 0;
    real sigma_adj = variance_adjusted_sgt(sigma, lambda, p, q);
    
    if (is_inf(q) && is_inf(p)) 
      return uniform_lpdf(x | mu - sigma_adj, mu + sigma_adj);
    
    real r = mean_centered_sgt(x, sigma_adj, lambda, p, q) - mu;
    real s = r < 0 ? -1 : 1;
    
    if (is_inf(q) && !is_inf(p)) {
      out = (abs(r) ./ (sigma_adj * (1 + lambda * s))) ^ p;
      return log(p) - log(2) - log(sigma_adj) - lgamma(1.0 / p) - out;
    } else {
      out = log1p(abs(r) ^ p ./ (q * sigma_adj ^ p * pow(1 + lambda * s, p)));
    }
    
    return N * (log(p) - log2() - log(sigma_adj) - log(q) / p - lbeta(1.0 / p, q)) - (1.0 / p + q) * out;
  }
  
  real skew_generalized_t_lcdf(real x, real mu, real sigma, real lambda, real p, real q) {
    if (sigma <= 0) 
      reject("sigma must be > 0 found sigma = ", sigma);
    
    if (lambda >= 1 || lambda <= -1) 
      reject("lambda must be between (-1, 1) found lambda = ", sigma);
    
    if (p <= 0) 
      reject("p must be > 0 found p = ", p);
    
    if (q <= 0) 
      reject("q must be > 0 found q = ", q);
    
    if (is_inf(x) && x < 0) 
      return 0;
    
    if (is_inf(x) && x > 0) 
      return 1;
    
    real sigma_adj = variance_adjusted_sgt(sigma, lambda, p, q);
    real x_cent = mean_centered_sgt(x, sigma_adj, lambda, p, q);
    real r = x_cent - mu;
    real lambda_new;
    real r_new;
    
    if (r > 0) {
      lambda_new = -lambda;
      r_new = -r;
    } else {
      lambda_new = lambda;
      r_new = r;
    }
    
    if (!is_inf(p) && is_inf(q) && !is_inf(x)) 
      return log_sum_exp([log1m(lambda_new) + log2(), log(lambda_new - 1) + log2() + beta_lcdf((r_new / (sigma * (1 + lambda_new))) ^ p | p, 1)]);
    if (is_inf(p) && !is_inf(x)) 
      return uniform_lcdf(x | mu, sigma);
    
    if (is_inf(x) && x < 0) 
      return 0;
    
    if (is_inf(x) && x > 0) 
      return 1;
    
    return log_sum_exp([log1m(lambda_new) + log2(),
                        log(lambda_new - 1) + log2() + beta_lcdf(1.0 / (1 + q * (sigma * (1 - lambda_new) / (-r_new)) ^ p) | 1.0 / p, q)]);
  }
'
