//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// target data is a matrix 'y' with 'K' columns and 'N' rows
// covariates are a matrix 'x' with 'J' columns and 'N' rows
data {
  int<lower=1> K;
  int<lower=1> J;
  int<lower=0> N;
  vector[J] x[N];
  vector[K] y[N];
}

// // The parameters accepted by the model.
// // KxJ matrix beta of regression parameters
// // covariance matrix Sigma
// parameters {
//   matrix[K, J] beta;
//   cov_matrix[K] Sigma;
// }
// 
// // The model to be estimated. We model the output
// // 'y' to be normally distributed with mean 'mu'
// // and standard deviation 'Sigma'.
// // where the mean mu is beta * x[n]
// model {
//   vector[K] mu[N];
//   for (n in 1:N)
//     mu[n] = beta * x[n];
//   y ~ multi_normal(mu, Sigma);
// }

// using a choleskz decomposition for the covariance matrix of y
parameters {
  matrix[K, J] beta;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
}

model {
  vector[K] mu[N];
  matrix[K, K] L_Sigma;

  for (n in 1:N)
    mu[n] = beta * x[n];

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  
  // priors for parameters
  to_vector(beta) ~ normal(0, 5);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);

  y ~ multi_normal_cholesky(mu, L_Sigma);
}
