# this is not working
# simulate data
K_Y = 2
L_Y = 4
N_Y = 5

K_X = 3
L_X = 5
N_X = 7

tau = 100
TT = 1000

set.seed(123)
# # get better prior for the diagonal elements in Sigma
Sigma_1 = diag(1, nrow = K_Y * K_X, ncol = K_Y * K_X)
mu_1 = rnorm(K_Y * K_X)
vecH_true = as.vector(matrix(1, nrow = K_Y, ncol = K_X))
mu_2 = rnorm(K_X^2)
Sigma_2 = diag(1, nrow = K_X^2, ncol = K_X^2)
vecB_true = as.vector(matrix(1, nrow = K_X, ncol = K_X))

# this is just a copy, not the actual values
vecH_true = rep(1, K_Y * K_X)
H_true = matrix(vecH_true, nrow = K_Y, ncol = K_X)
vecB_true = diag(1, K_X)
B_true = matrix(vecB_true, nrow = K_X, ncol = K_X)

C_true = 1 * diag(K_Y)
D_true = 1 * diag(K_X)

b_true = matrix(NA, nrow = K_X, ncol = 1000)
b_true[,1] = 0
for (i in 2:1000) {
  b_true[,i] = mvnfast::rmvn(1, B_true %*% b_true[,i - 1], D_true)
}

V_true = prcomp(matrix(rnorm(N_Y * TT), nrow = TT))$rotation
V_true = diag(N_Y:1)
X = V_true[,1:K_X] %*% b_true
X = X + matrix(rnorm(length(X), sd = 0.1), nrow = nrow(X), ncol = ncol(X))

a_true = matrix(NA, nrow = K_Y, ncol = 1000)
for (i in 1:1000) {
  a_true[,i] = mvnfast::rmvn(1, H_true %*% b_true[,i], C_true)
}

eigen_sie_detrended = prcomp(sie_detrended[complete.cases(sie_detrended),], scale. = FALSE, center = TRUE)
U_true = eigen_sie_detrended$rotation
a_true = eigen_sie_detrended$x[, 1:K_Y] %>% t()
Y = U_true[,1:K_Y] %*% a_true
Y = Y + matrix(rnorm(length(Y), sd = 0.1), nrow = nrow(Y), ncol = ncol(Y))

true_Y = Y
# set first hundred values to be missing
Y[,1:99] = NA


# do an eigenvector decomposition
eigen_Y = prcomp(t(Y)[complete.cases(t(Y)),], scale. = FALSE, center = FALSE)
# use the first two principal components, see
plot(eigen_Y$sdev)
eigen_Y$sdev / sum(eigen_Y$sdev)
K_Y = 2
L_Y = 4
# get dimensionality of Y
N_Y = length(eigen_Y$sdev)

# get basis vectors in matrix U
U = eigen_Y$rotation[,1:K_Y]
tildeU = eigen_Y$rotation[,(K_Y + 1):L_Y]


# total sample variance after accounting for the first K_Y EOFs
L_Y = 4
c = sum(eigen_Y$sdev[(L_Y + 1):N_Y])
tildeLambda = diag(-(eigen_Y$sdev[(K_Y + 1):L_Y]) / (c * (c + eigen_Y$sdev[(K_Y + 1):L_Y])))


# do an eigenvector decomposition on covariates
eigen_X = prcomp(t(X))
# use the first 3 principal components, see
plot(eigen_X$sdev)
eigen_X$sdev / sum(eigen_X$sdev)
K_X = 3
N_X = length(eigen_X$sdev)
L_X = 5

# get basis vectors in matrix U
V = eigen_X$rotation[,1:K_X]
tildeV = eigen_X$rotation[,(K_X + 1):L_X]
# get amplitudes in matrix b
b = t(V) %*% X
b[,1:2]

tildeLambda_X = diag(-(eigen_X$sdev[(K_X + 1):L_X]) / (c * (c + eigen_X$sdev[(K_X + 1):L_X])))

# these do not hold any sway
r_true = 1
s_true = 1
tildeR_true = diag(c, N_Y) + tildeU %*% diag(eigen_Y$sdev[(K_Y + 1):L_Y]) %*% t(tildeU)
R_true = r_true * tildeR_true
tildeS_true = diag(c, N_X) + tildeV %*% diag(eigen_X$sdev[(K_X + 1):L_X]) %*% t(tildeV)
S_true = s_true * tildeS_true


# first observed value in Y
tau = 100
# total number of observations
TT = ncol(Y)

################## Gibbs sampler ################## 
# 1000 iterations take about 5 minutes
niter = 100
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
{
  # seems reasonable
  # a.prior = t(U) %*% Y
  # b.prior = t(V) %*% X
  a.prior = matrix(rnorm(TT * K_Y), nrow = K_Y)
  b.prior = matrix(rnorm(TT * K_X), nrow = K_X)
  a = a.prior
  b = b.prior
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
}

##### Gibbs Sampler ####
# compute constants
{
  outside_helper = t(U) %*% Y
  outside_helper2 = t(V) %*% X
  
  E_Y = diag(1 / c, N_Y) - tildeU %*% tildeLambda %*% t(tildeU)
  E_X = diag(1 / c, N_X) - tildeV %*% tildeLambda_X %*% t(tildeV)
}
# profvis::profvis({
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
  for (t in tau:TT) {
    a_new[,t] = mvnfast::rmvn(n = 1, mu = a.mu[,t], sigma = a.Sigma)
  }
  
  for (t in 1:(tau - 1)) {
    a_new[,t] = mvnfast::rmvn(n = 1, mu = H %*% b[,t], sigma = C)
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
    b_new[,t] = mvnfast::rmvn(n = 1, mu = b.mu[,t], sigma = b.Sigma)
  }
  
  # compute mu and Sigma for last b_T
  b.SigmaT = solve(1 / (s * c) * diag(K_X) +
                     t(H) %*% solve(C) %*% H + 
                     solve(D))
  b.mu[,TT] = b.SigmaT %*% (1 / (s * c) * outside_helper2[,TT] + 
                              t(H) %*% solve(C) %*% a[,t] + 
                              solve(D) %*% B %*% b[, t - 1])
  b_new[,TT] = mvnfast::rmvn(n = 1, mu = b.mu[,TT], sigma = b.SigmaT)
  
  # compute mu and Sigma for b_1
  # TODO get better priors for b_1
  mu.b = rnorm(K_X, 0, 1)
  Sigma_b = diag(K_X)
  b.Sigma1 = solve(t(B) %*% solve(D) %*% B + solve(Sigma_b))
  b.mu[,1] = b.Sigma1 %*% (t(B) %*% solve(D) %*% b.mu[,2] + solve(Sigma_b) %*% mu.b)
  b_new[,1] = mvnfast::rmvn(n = 1, mu = b.mu[,1], sigma = b.Sigma1)
  # add to results
  b = b_new
  results$b[[iter]] = b_new
  
  
  ################################################################################
  # # transition matrix H
  # compute the inverse first, as this is easier to sum up
  vecH.Sigma_inv = solve(Sigma_1)
  # compute a helper for mu
  vecH.mu_helper = solve(Sigma_1) %*% mu_1
  for (t in tau:TT) {
    J = kronecker(t(b[,t]), diag(1, K_Y))
    JtCinv = t(J) %*% solve(C)
    vecH.Sigma_inv = vecH.Sigma_inv + JtCinv %*% J
    vecH.mu_helper = vecH.mu_helper + JtCinv %*% a[,t]
    # print(t)
    # print(vecH.mu_helper)
  }
  # then invert it to get actual covariance matrix
  vecH.Sigma = solve(vecH.Sigma_inv)
  vecH.mu = vecH.Sigma %*% vecH.mu_helper
  
  # sample actual vecH and H
  vecH = mvnfast::rmvn(n = 1, mu = vecH.mu, sigma = vecH.Sigma)[1,]
  H = matrix(vecH, nrow = K_Y, ncol = K_X)
  # add to results
  results$vecH[[iter]] = vecH
  results$H[[iter]] = H
  
  ################################################################################
  # # transition matrix B
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
  vecB = mvnfast::rmvn(n = 1, mu = vecB.mu, sigma = vecB.Sigma)[1,]
  B = matrix(vecB, nrow = K_X, ncol = K_X)
  
  # add to results
  results$vecB[[iter]] = vecB
  results$B[[iter]] = B
  
  ################################################################################
  # # covariances of data model r and s
  
  # TODO what are T_X and T_Y?
  T_X = TT
  T_Y = TT - tau
  
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
# })

str(results)

pl = c(50, 75, 100)

get_RMSE = function(results, true_Y, iter) {
  true_Y = t(true_Y)
  res_data = data.frame(iter = iter, 
                        King_Hakon = NA, Ross = NA, East_Antarctica = NA, 
                        Weddell = NA, Bellingshausen_Amundsen_Sea = NA)
  for (i in iter) {
    Y_est = (results$U %*% results$a[[i]]) %>% t() %>% as.data.frame()
    res_data[i, 2:6] = sqrt(colMeans(((Y_est - true_Y)^2)[1:tau,]))
  }
  res_data
}
rmse_data = get_RMSE(results, true_Y, results$iter)
plot(rmse_data$King_Hakon, type = "l")
plot(rmse_data$Ross, type = "l")
plot(rmse_data$East_Antarctica, type = "l")
plot(rmse_data$Weddell, type = "l")
plot(rmse_data$Bellingshausen_Amundsen_Sea, type = "l")

# check results graphically
plot(results$a[[pl[1]]][1,50:(2 * tau)], type = "l", col = 4)
lines(results$a[[pl[2]]][1,50:(2 * tau)], col = 2)
lines(results$a[[pl[3]]][1,50:(2 * tau)], col = 3)
lines(-a_true[1,50:(2 * tau)])

cors = lapply(results$a, function(x) {
  cor(x[1,1:tau], a_true[1,1:tau])
}) %>%
  unlist()
cors
hist(cors)

plot(results$b[[pl[1]]][1,], type = "l", col = 4)
lines(results$b[[pl[2]]][1,], col = 2)
lines(results$b[[pl[3]]][1,], col = 3)
lines(b_true[1,])

# vecH is very weak looking
plot(results$vecH[[pl[1]]], type = "l", col = 4)
lines(results$vecH[[pl[2]]], col = 2)
lines(results$vecH[[pl[3]]], col = 3)
lines(vecH_true)

# par(mfrow = c(1, 3))
# image(results$H[[pl[1]]])
# image(results$H[[pl[2]]])
# image(results$H[[pl[3]]])
# par(mfrow = c(1, 1))


plot(results$vecB[[pl[1]]], type = "l", col = 4)
lines(results$vecB[[pl[2]]], col = 2)
lines(results$vecB[[pl[3]]], col = 3)
lines(vecB_true)

# par(mfrow = c(1, 3))
# image(results$B[[pl[1]]])
# image(results$B[[pl[2]]])
# image(results$B[[pl[3]]])
# par(mfrow = c(1, 1))

# r is not really that close
plot(results$r[-(1:10)], type = "l")
abline(h = r_true)
# s getting closer
plot(results$s[-(1:10)], type = "l")
abline(h = s_true)

plot(as.vector(results$C[[pl[1]]]), type = "l", col = 4)
lines(as.vector(results$C[[pl[2]]]), col = 2)
lines(as.vector(results$C[[pl[3]]]), col = 3)
lines(C_true)

plot(as.vector(results$D[[pl[1]]]), type = "l", col = 4)
lines(as.vector(results$D[[pl[2]]]), col = 2)
lines(as.vector(results$D[[pl[3]]]), col = 3)
lines(D_true)

par(mfrow = c(1, 4))
image(results$D[[pl[1]]])
image(results$D[[pl[2]]])
image(results$D[[pl[3]]])
image(D_true)
par(mfrow = c(1, 1))
