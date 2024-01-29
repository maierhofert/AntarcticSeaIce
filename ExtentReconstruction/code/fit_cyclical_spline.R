library("splines")
library("rstan")
set.seed(1234)
X <- seq(from=0, to=365, by=1) # generating inputs
# B <- t(bs(X, knots=seq(366 * 1/4,366 * 3/4,length.out = 2), Boundary.knots = c(0, 366), degree=3, intercept = TRUE)) # creating the B-splines
# str(B)
# dim(B)

B = t(mgcv::cSplineDes(x = X, knots = seq(0,366, length.out = 5), ord = 3))
dim(B)

num_data <- length(X); 
num_basis <- nrow(B)

a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
Y_true <- as.vector(a%*%B) # generating the output
Y <- Y_true + rnorm(length(X),0,.2) # adding noise
rstan_options(auto_write = TRUE)

plot(X, Y_true, type = "l")
points(X, Y)
mean(Y_true)

# options(mc.cores = parallel::detectCores())
sm<-stan_model("fit_cyclical_spline.stan")
fit<-sampling(sm,iter=500, 
              data = list(X = X, B = B, num_data = num_data, num_basis = num_basis, 
                          a = a, Y_true = Y_true, Y = Y),
              control = list(adapt_delta=0.95))

plot(X, Y_true, type = "l")
points(X, Y)

# estimates
fit_summary = summary(fit, pars = "Y_hat")
fit_summary$summary[,"mean"]
lines(X, fit_summary$summary[,"mean"], col = 2)

