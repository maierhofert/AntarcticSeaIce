library(rstan)

options(mc.cores = 10)
rstan_options(auto_write = TRUE)


stan_data = list(K = 5,
                 J = 10,
                 N = 501,
                 x = all_dat[11000:11500,200:209],
                 y = all_dat[11000:11500,4:8])

multivariate_regression = stan_model("multivariate_regression.stan")
# horseshoe prior can be implemented like this:
# https://discourse.mc-stan.org/t/horseshoe-prior-in-rstan/14561
# or regularized horseshoe here?
# https://discourse.mc-stan.org/t/bayes-sparse-regression-reg-horseshoe-for-multinomial-model/4336/3

# cyclical spline functions like this

# actually fit the model
# this is very small in order to keep it fast
fit1 <- sampling(multivariate_regression, data = stan_data, 
                 chains = 2, iter = 100)

print(fit1, pars = "beta", probs = c(0.025, 0.5, 0.975))
traceplot(fit1)
plot(fit1)

