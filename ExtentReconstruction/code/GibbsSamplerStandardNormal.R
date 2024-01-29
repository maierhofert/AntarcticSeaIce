# make up target variable
day = 1:100
set.seed(123)
sie = matrix(rnorm(5 * length(day)), ncol = 5)
sie[1:10,] = NA
# covariate
# make up covariates
cov = matrix(rnorm(6 * length(day)), ncol = 6)

# define target variable and covariate
Y = sie
X = cov

# impute missing values
# for now just use the previous observation for simplicity
# Y[359,] = Y[358,]
# standardize Y
Y = scale(Y, center = TRUE, scale = FALSE)

# do an eigenvector decomposition
eigen_Y = prcomp(Y[complete.cases(Y),], scale. = FALSE, center = FALSE)
# use the first two principal components, see
plot(eigen_Y$sdev)
K_Y = 2
L_Y = 4
# get dimensionality of Y
N_Y = ncol(Y)

# get basis vectors in matrix U
U = eigen_Y$rotation[,1:K_Y]
tildeU = eigen_Y$rotation[,(K_Y + 1):L_Y]

# transpose Y to make it N_Y by T
Y = t(Y)
# get amplitudes in matrix a
a1 = Y[,1] %*% U
a = t(U) %*% Y

dim(U %*% a)
eta = U %*% a - Y
summary(eta)
apply(eta, 1, sd)
# total sample variance after accounting for the first K_Y EOFs
L_Y = 4
c = sum(eigen_Y$sdev[(L_Y + 1):N_Y])
tildeLambda = diag(-(eigen_Y$sdev[(K_Y + 1):L_Y]) / (c * (c + eigen_Y$sdev[(K_Y + 1):L_Y])))


# # create scaled values
# Y %*% eigen_Y$rotation - eigen_Y$x
# # recreate original
# Y %*% eigen_Y$rotation %*% t(eigen_Y$rotation) - Y
# eigen_Y$x %*% t(eigen_Y$rotation) - Y



# do an eigenvector decomposition on covariates
eigen_X = prcomp(X)
# use the first 3 principal components, see
plot(eigen_X$sdev)
K_X = 3
N_X = length(eigen_X$sdev)
L_X = 5

# get basis vectors in matrix U
V = eigen_X$rotation[,1:K_X]
tildeV = eigen_X$rotation[,(K_X + 1):L_X]
# get amplitudes in matrix b
X = t(X)
b = t(V) %*% X
b[,1:2]

tildeLambda_X = diag(-(eigen_X$sdev[(K_X + 1):L_X]) / (c * (c + eigen_X$sdev[(K_X + 1):L_X])))


# this is 93 in original paper
tau = 11
TT = ncol(Y)

################## Gibbs sampler ################## 
niter = 3
iter = 1
results = list(iter = 1:niter,
               Y = Y,
               U = U, tildeU = tildeU,
               N_Y = nrow(Y), K_Y = dim(U)[2], 
               L_Y = 4, c = c,
               X = X,
               V = V, tildeV = tildeV,
               N_X = N_X, K_X = K_X,
               tau = 11, TT = ncol(Y),
               a = list(),
               b = list(),
               vecH = list(), H = list(),
               vecB = list(), B = list(),
               r = rep(NA, niter), tildeR = NA,
               s = rep(NA, niter), tildeS = NA,
               C = list(), D = list())

##### priors ####

# seems reasonable
a.prior = t(U) %*% Y
b.prior = t(V) %*% X
# a.priori = matrix(rnorm(TT * K_Y), ncol = K_Y)
# b.priori = matrix(rnorm(TT * K_X), ncol = K_X)
results$a[[1]] = a.prior
results$b[[1]] = b.prior

# # get better prior for the diagonal elements in Sigma
Sigma_1 = diag(1, nrow = K_Y * K_X, ncol = K_Y * K_X)
mu_1 = rnorm(K_Y * K_X)
vecH = MASS::mvrnorm(n = 1, 
                     mu = mu_1, 
                     Sigma = Sigma_1)
mu_2 = rnorm(K_X^2)
Sigma_2 = diag(1, nrow = K_X^2, ncol = K_X^2)
vecB = MASS::mvrnorm(n = 1, 
                     mu = mu_2, 
                     Sigma = Sigma_2)
# this is just a copy, not the actual values
H = matrix(vecH, nrow = K_Y, ncol = K_X)
B = matrix(vecB, nrow = K_X, ncol = K_X)
# add to results
results$vecH[[1]] = vecH
results$vecB[[1]] = vecB
results$H[[1]] = H
results$B[[1]] = B

# good priors
r = invgamma::rinvgamma(1, 2)
s = invgamma::rinvgamma(1, 2)
tildeR = diag(c, N_Y) + tildeU %*% diag(eigen_Y$sdev[(K_Y + 1):L_Y]) %*% t(tildeU)
R = r * tildeR
tildeS = diag(c, N_X) + tildeV %*% diag(eigen_X$sdev[(K_X + 1):L_X]) %*% t(tildeV)
S = s * tildeS
# add to results
results$r[1] = r
results$s[1] = s
results$tildeR = tildeR
results$tildeS = tildeS

# get better value for prefactor of Sigma
C = CholWishart::rInvWishart(n = 1L, df = K_Y + 2, Sigma = 1 * diag(K_Y))[,,1]
D = CholWishart::rInvWishart(n = 1L, df = K_X + 2, Sigma = 1 * diag(K_X))[,,1]
# add to results
results$C[[1]] = C
results$D[[1]] = D

##### Gibbs Sampler ####

outside_helper = t(U) %*% Y
outside_helper2 = t(V) %*% X

for (iter in 2:niter) {
  ################################################################################
  # target amplitudes a_t
  a.Sigma = solve(1 / (r * c) * diag(K_Y) + solve(C))
  a.mu = matrix(NA, ncol = TT, nrow = K_Y)
  for (t in tau:TT) {
    a.mu[,t] = a.Sigma %*% (1 / (r * c) * outside_helper[,t] + solve(C) %*% H %*% b[,t])
  }
  # sample actual a_t
  a_new = matrix(NA, nrow = K_Y, ncol = TT)
  for (t in 1:TT) {
    a_new[,t] = MASS::mvrnorm(n = 1, mu = a.mu[,t], Sigma = a.Sigma)
  }
  for (t in 1:(tau - 1)) {
    a_new[,t] = MASS::mvrnorm(n = 1, mu = H %*% b[,t], Sigma = C)
  }
  # add to results
  a = a_new
  results$a[[iter]] = a_new
  # plot(a_new[1,])
  # points(a[1,], col = "red")
  # points(a_new[1,] - a[1,], col = "blue")
  
  ################################################################################
  # ### covariate amplitudes b_t #####################
  
  # compute mu and Sigma for b_2 to b_{TT-1}
  b.Sigma = solve(1 / (s * c) * diag(K_X) +
                    t(H) %*% solve(C) %*% H + 
                    solve(D) + 
                    t(B) %*% solve(D) %*% B)
  
  b.mu = matrix(NA, ncol = TT, nrow = K_X)
  for (t in 2:(TT - 1)) {
    b.mu[,t] = b.Sigma %*% (1 / (s * c) * outside_helper2[,t] + 
                              t(H) %*% solve(C) %*% a[,t] + 
                              solve(D) %*% B %*% b[, t - 1] +
                              t(B) %*% solve(D) %*% b[, t + 1])
  }
  # sample actual b_t
  b_new = matrix(NA, nrow = K_X, ncol = TT)
  for (t in 2:(TT - 1)) {
    b_new[,t] = MASS::mvrnorm(n = 1, mu = b.mu[,t], Sigma = b.Sigma)
  }
  
  # compute mu and Sigma for last b_T
  b.SigmaT = solve(1 / (s * c) * diag(K_X) +
                     t(H) %*% solve(C) %*% H + 
                     solve(D))
  b.mu[,TT] = b.SigmaT %*% (1 / (s * c) * outside_helper2[,TT] + 
                              t(H) %*% solve(C) %*% a[,t] + 
                              solve(D) %*% B %*% b[, t - 1])
  b_new[,TT] = MASS::mvrnorm(n = 1, mu = b.mu[,TT], Sigma = b.SigmaT)
  
  # compute mu and Sigma for b_1
  # TODO get better priors for b_1
  mu.b = rnorm(K_X, 0, 1)
  Sigma_b = diag(K_X)
  b.Sigma1 = solve(t(B) %*% solve(D) %*% B + solve(Sigma_b))
  b.mu[,1] = b.Sigma1 %*% (t(B) %*% solve(D) %*% b.mu[,2] + solve(Sigma_b) %*% mu.b)
  b_new[,1] = MASS::mvrnorm(n = 1, mu = b.mu[,1], Sigma = b.Sigma1)
  # add to results
  b = b_new
  results$b[[iter]] = b_new
  
  
  ################################################################################
  # # transition matrix H
  # TODO figure out when tau actually is
  tau = 11
  
  # compute the inverse first, as this is easier to sum up
  vecH.Sigma_inv = solve(Sigma_1)
  # compute a helper for mu
  vecH.mu_helper = solve(Sigma_1) %*% mu_1
  for (t in tau:TT) {
    J = kronecker(t(b[,t]), diag(1, K_Y))
    JtCinv = t(J) %*% solve(C)
    vecH.Sigma_inv = vecH.Sigma_inv + JtCinv %*% J
    vecH.mu_helper = vecH.mu_helper + JtCinv %*% a[,t]
  }
  # then invert it to get actual covariance matrix
  vecH.Sigma = solve(vecH.Sigma_inv)
  vecH.mu = vecH.Sigma %*% vecH.mu_helper
  
  # sample actual vecH and H
  vecH = MASS::mvrnorm(n = 1, mu = vecH.mu, Sigma = vecH.Sigma)
  H = matrix(vecH, nrow = K_Y, ncol = K_X)
  # add to results
  results$vecH[[iter]] = vecH
  results$H[[iter]] = H
  
  ################################################################################
  # # transition matrix B
  # TODO figure out when tau actually is
  
  # compute the inverse first, as this is easier to sum up
  vecB.Sigma_inv = solve(Sigma_2)
  # compute a helper for mu
  vecB.mu_helper = solve(Sigma_2) %*% mu_2
  for (t in tau:TT) {
    tildeJ_t_minus_1 = kronecker(t(b[,t - 1]), diag(1, K_X))
    JtDinv = t(tildeJ_t_minus_1) %*% solve(D)
    vecB.Sigma_inv = vecB.Sigma_inv + JtDinv %*% tildeJ_t_minus_1
    vecB.mu_helper = vecB.mu_helper + JtDinv %*% b[,t]
  }
  # then invert it to get actual covariance matrix
  vecB.Sigma = solve(vecB.Sigma_inv)
  vecB.mu = vecB.Sigma %*% vecB.mu_helper
  
  # sample actual vecH and H
  vecB = MASS::mvrnorm(n = 1, mu = vecB.mu, Sigma = vecB.Sigma)
  B = matrix(vecB, nrow = K_X, ncol = K_X)
  
  # add to results
  results$vecB[[iter]] = vecB
  results$B[[iter]] = B
  
  ################################################################################
  # # covariances of data model r and s
  E_Y = diag(1 / c, N_Y) - tildeU %*% tildeLambda %*% t(tildeU)
  E_X = diag(1 / c, N_X) - tildeV %*% tildeLambda_X %*% t(tildeV)
  
  # TODO what are T_X and T_Y?
  T_X = TT
  T_Y = TT
  
  beta_1_new = 2
  for (t in tau:TT) {
    beta_1_new = beta_1_new + 1 / 2 * t(Y[,t] - U %*% a[,t]) %*% E_Y %*% (Y[,t] - U %*% a[,t])
  }
  r_new = invgamma::rinvgamma(n = 1, T_Y * N_Y / 2 + 1, beta_1_new)
  
  beta_2_new = 2
  for (t in tau:TT) {
    beta_2_new = beta_2_new + 1 / 2 * t(X[,t] - V %*% b[,t]) %*% E_X %*% (X[,t] - V %*% b[,t])
  }
  s_new = invgamma::rinvgamma(n = 1, T_X * N_X / 2 + 1, beta_2_new)
  
  # add to results
  r = r_new
  results$r[iter] = r_new
  s = s_new
  results$s[iter] = s_new
  
  ################################################################################
  # # covariance matrices of data model C and D
  # TODO add better diagonal value
  W_C = diag(1, 1 * K_Y)
  W_D = diag(1, 1 * K_X)
  
  nu_Y = K_Y + 2
  nu_X = K_X + 2
  
  C_mean = W_C
  for (t in tau:TT) {
    C_mean = C_mean + (a[,t] - H %*% b[,t]) %*% t(a[,t] - H %*% b[,t])
  }
  C_new = CholWishart::rInvWishart(n = 1L, 
                                   df = nu_Y + TT - tau + 1, 
                                   Sigma = C_mean)[,,1]
  
  
  D_mean = W_D
  for (t in tau:TT) {
    D_mean = D_mean + (b[,t] - B %*% b[,t - 1]) %*% t(b[,t] - B %*% b[,t - 1])
  }
  D_new = CholWishart::rInvWishart(n = 1L, 
                                   df = nu_X + TT - tau + 1, 
                                   Sigma = D_mean)[,,1]
  
  # add to results
  C = C_new
  results$C[[iter]] = C_new
  D = D_new
  results$D[[iter]] = D_new
  
}

str(results)

# check results graphically
plot(results$a[[1]][1,], type = "l")
lines(results$a[[2]][1,], col = 2)
lines(results$a[[3]][1,], col = 3)

plot(results$vecH[[1]], type = "l")
lines(results$vecH[[2]], col = 2)
lines(results$vecH[[3]], col = 3)

plot(results$vecB[[1]], type = "l")
lines(results$vecB[[2]], col = 2)
lines(results$vecB[[3]], col = 3)

plot(results$r, type = "l")
plot(results$s, type = "l")

plot(as.vector(results$C[[1]]), ylim = c(0, 20))
points(as.vector(results$C[[2]]), col = 2)
points(as.vector(results$C[[3]]), col = 3)

plot(as.vector(results$D[[1]]), ylim = c(0, 3))
points(as.vector(results$D[[2]]), col = 2)
points(as.vector(results$D[[3]]), col = 3)


