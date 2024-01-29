# read Mark's data

# TODO read in required covariates
# read in covariates
amelia_ghcnd_groupASN1 = readRDS("../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN1.RDS")
amelia_ghcnd_groupASN2 = readRDS("../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN2.RDS")
amelia_ghcnd_groupASN3 = readRDS("../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN3.RDS")
amelia_ghcnd_groupASN4 = readRDS("../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN4.RDS")

cov1 = amelia_ghcnd_groupASN1$imputations[["imp1"]][,c(22, 4:21)]
cov2 = amelia_ghcnd_groupASN2$imputations[["imp1"]][,4:37]
cov3 = amelia_ghcnd_groupASN3$imputations[["imp1"]][,4:21]
cov4 = amelia_ghcnd_groupASN4$imputations[["imp1"]][,4:22]
cov = cbind(cov1, cov2, cov3, cov4)
for (i in 1:ncol(cov)) {
  cov[,i] = c(cov[,i])
}
summary(cov)

# smooth the covariates before putting them into the model
i = 2
cov_s = cov
for (i in 2:ncol(cov)) {
  cov_s[,i] = ksmooth(cov$Date, cov[,i], kernel = "normal", bandwidth = 30.4)$y
}
plot(cov$Date, cov[,2], xlim = c(lubridate::ymd("1901-01-01"), lubridate::ymd("1910-01-01")), type = "l")
lines(cov$Date, cov_s[,2], col = 2)
# is this a problem? not unit variance anymore
apply(cov_s, 2, sd)

# add lagged versions of all covariates
lags = c(2, 7, 29, 61)
col_names = names(cov)
jmax = ncol(cov)
for (j in 2:jmax) {
  for (i in seq_along(lags)) {
    cov[,paste0(col_names[j], "_l", lags[i])] = lag(cov[,col_names[j]], n = lags[i])
  }
}

all_data = full_join(cbind(data.frame(tdate, Doy, Date = sie_date),
                           sie_detrended),
                     cov, 
                     by = "Date")
all_data = all_data[order(all_data$Date),]

min(which(complete.cases(all_data$King_Hakon)))

# define target variable and covariate
Y = all_data[34500:43464, 
             c("King_Hakon", "Ross", "East_Antarctica", 
               "Weddell", "Bellingshausen_Amundsen_Sea")]
which(!complete.cases(Y))
# there is a few days missing, fill that in for now
Y[which(!complete.cases(Y)), ] = Y[which(!complete.cases(Y)) - 1, ]
true_Y = Y
# set first hundred values to be missing
Y[1:365,] = NA
# this is a small training data set
X = all_data[34500:43464, 9:453]


# do an eigenvector decomposition
eigen_Y = prcomp(Y[complete.cases(Y),], scale. = FALSE, center = FALSE)
K_Y = 3
L_Y = 4
# get dimensionality of Y
N_Y = ncol(Y)

# use the first two principal components, see
plot(eigen_Y$sdev)
eigen_Y$sdev / sum(eigen_Y$sdev)
sum(eigen_Y$sdev[1:K_Y]) / sum(eigen_Y$sdev)


# get basis vectors in matrix U
U = eigen_Y$rotation[,1:K_Y]
tildeU = eigen_Y$rotation[,(K_Y + 1):L_Y]

# transpose Y to make it N_Y by T
Y = t(Y)
# # get amplitudes in matrix a
# a1 = Y[,1] %*% U
a = t(U) %*% Y

# dim(U %*% a)
# eta = U %*% a - Y
# summary(eta)
# apply(eta, 1, sd)
# total sample variance after accounting for the first K_Y EOFs
L_Y = 4
c = sum(eigen_Y$sdev[(L_Y + 1):N_Y])
tildeLambda = diag(-(eigen_Y$sdev[(K_Y + 1):L_Y]) / (c * (c + eigen_Y$sdev[(K_Y + 1):L_Y])),
                   nrow = L_Y - K_Y, ncol = L_Y - K_Y)


# # create scaled values
# Y %*% eigen_Y$rotation - eigen_Y$x
# # recreate original
# Y %*% eigen_Y$rotation %*% t(eigen_Y$rotation) - Y
# eigen_Y$x %*% t(eigen_Y$rotation) - Y



# do an eigenvector decomposition on covariates
# maybe do a smoothed version of this
# eigen_X = prcomp(cov[-(1:61),-1])
eigen_X = prcomp(X)
# use the first ??? principal components, see
plot(eigen_X$sdev)
abline(h = 1)
sum(eigen_X$sdev > 1)

image(cor(X))
image(abs(cor(X)))

# there is still structure in the first 20 or so I think
image(abs(eigen_X$rotation[,1:10]), xlab = "Covariates", ylab = "Eigenvectors", zlim = c(0, 0.3))
image(abs(eigen_X$rotation[,c(11:20)]), xlab = "Covariates", ylab = "Eigenvectors", zlim = c(0, 0.3))
image(abs(eigen_X$rotation[,21:30]), xlab = "Covariates", ylab = "Eigenvectors", zlim = c(0, 0.3))
# but just white noise after that
image(abs(eigen_X$rotation[,31:39]), xlab = "Covariates", ylab = "Eigenvectors", zlim = c(0, 0.3))

# There are sets of two with opposing loadings, I assume that's TMIN and TMAX
# this reappears in stronger effects for every other PC in regression
image(abs(eigen_X$rotation[,1:25]), xlab = "Covariates", ylab = "Eigenvectors", zlim = c(0, 0.3))

image(abs(eigen_X$rotation[,seq(1,80, by = 5)]), xlab = "Covariates", ylab = "Eigenvectors", zlim = c(0, 0.3))


K_X = 22
abline(v = K_X)
sum(eigen_X$sdev[1:K_X]) / sum(eigen_X$sdev)

N_X = length(eigen_X$sdev)
L_X = 25

# get basis vectors in matrix U
V = eigen_X$rotation[,1:K_X]
tildeV = eigen_X$rotation[,(K_X + 1):L_X]
# get amplitudes in matrix b
X = t(X)
b = t(V) %*% X
b[,1:2]


tildeLambda_X = diag(-(eigen_X$sdev[(K_X + 1):L_X]) / (c * (c + eigen_X$sdev[(K_X + 1):L_X])),
                     nrow = L_X - K_X, ncol = L_X - K_X)


# first observed value in Y
tau = 366
Y_miss = !complete.cases(t(Y))
Y_miss_total = sum(Y_miss)
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
  
  tildeR = diag(c, N_Y) + tildeU %*% 
    diag(eigen_Y$sdev[(K_Y + 1):L_Y], nrow = L_Y - K_Y, ncol = L_Y - K_Y) %*% 
    t(tildeU)
  R = r * tildeR
  tildeS = diag(c, N_X) + tildeV %*% 
    diag(eigen_X$sdev[(K_X + 1):L_X], nrow = L_X - K_X, ncol = L_X - K_X) %*% 
    t(tildeV)
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
  # if (iter %% 10 == 0) print(iter)
  ################################################################################
  # target amplitudes a_t
  a.Sigma = solve(1 / (r * c) * diag(K_Y) + solve(C))
  a.mu = matrix(NA, ncol = TT, nrow = K_Y)
  # for (t in which(!Y_miss)) {
  for (t in tau:TT) {
    a.mu[,t] = a.Sigma %*% (1 / (r * c) * outside_helper[,t] + solve(C) %*% H %*% b[,t])
  }
  # sample actual a_t
  a_new = matrix(NA, nrow = K_Y, ncol = TT)
  for (t in tau:TT) {
    # TODO break out the covariance matrix and multiply it on afterwards
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
    # # this is not actually faster
    # b_new[,t] = mvnfast::rmvn(n = 1, mu = b.mu[,t], sigma = diag(K_X))
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
# # save result
# saveRDS(results, "../data/GibbsSamplerResults/results.RDS")
# results = readRDS("../data/GibbsSamplerResults/results.RDS")

pl = c(50, 75, 100)

# compute estimated Y
results$Y_mean = list()
results$Y_draw = list()
for (iter in results$iter) {
  Y_est = (results$U %*% results$a[[iter]]) %>% t()
  results$Y_est[[iter]] = Y_est
  Y_draw = apply(Y_est, 1, function(mu) {
    mvnfast::rmvn(n = 1, mu = mu, sigma = results$r[i] * results$tildeR)
  }) %>% t()
  results$Y_draw[[iter]] = Y_draw
}
# plot last draw from sampler
plot(Y_draw[1:(2*tau),1], type = "l", col = 2, main = "King Hakon")
lines(Y_est[1:(2*tau),1], col = 4)
lines(true_Y[1:(2*tau),1], col = 1)

plot(Y_draw[1:(2*tau),2], type = "l", col = 2, main = "Ross")
lines(Y_est[1:(2*tau),2], col = 4)
lines(true_Y[1:(2*tau),2], col = 1)

plot(Y_draw[1:(2*tau),3], type = "l", col = 2, main = "East_Antarctica")
lines(Y_est[1:(2*tau),3], col = 4)
lines(true_Y[1:(2*tau),3], col = 1)

plot(Y_draw[1:(2*tau),4], type = "l", col = 2, main = "Weddell")
lines(Y_est[1:(2*tau),4], col = 4)
lines(true_Y[1:(2*tau),4], col = 1)

plot(Y_draw[1:(2*tau),5], type = "l", col = 2, main = "Bellingshausen-Amundsen")
lines(Y_est[1:(2*tau),5], col = 4)
lines(true_Y[1:(2*tau),5], col = 1)


# compute correlation between true and estimated Y
correlations = lapply(results$Y_est[-(1:10)], function(x) {
  diag(cor(x[1:tau,], true_Y[1:tau, ]))
}) %>% unlist() %>% matrix(ncol = 5)
summary(correlations)
hist(correlations[,1])
hist(correlations[,2])
hist(correlations[,3])
hist(correlations[,4])
hist(correlations[,5])


# compute RMSE
get_RMSE = function(results, true_Y, iter) {
  res_data = data.frame(iter = iter, 
                        King_Hakon = NA, Ross = NA, East_Antarctica = NA, 
                        Weddell = NA, Bellingshausen_Amundsen_Sea = NA)
  
  for (i in seq_along(iter)) {
    # Y_est = (results$U %*% results$a[[iter[i]]]) %>% t() %>% as.data.frame()
    Y_est = results$Y_est[[iter[i]]]
    # Y_draw = results$Y_draw[[iter[i]]]
    res_data[i + 1, 2:6] = sqrt(colMeans(((Y_est - true_Y)^2)[1:tau,]))
  }
  res_data
}
rmse_data = get_RMSE(results, true_Y, results$iter[-(1:10)])
# compare to estimating true data mean or sample mean 
rmse_base2 = sqrt(colMeans(((0 - true_Y)^2)[1:tau,]))
rmse_base = sqrt(colMeans(((rowMeans(Y, na.rm = TRUE) - true_Y)^2)[1:tau,]))

# summary(Y_est)
# summary(true_Y)

library(ggplot2)
ggplot(rmse_data) +
  geom_line(aes(x = iter, y = King_Hakon)) +
  geom_hline(yintercept = rmse_base["King_Hakon"], col = "red", lty = "dashed")
ggplot(rmse_data) +
  geom_line(aes(x = iter, y = Ross)) +
  geom_hline(yintercept = rmse_base["Ross"], col = "red", lty = "dashed")
ggplot(rmse_data) +
  geom_line(aes(x = iter, y = East_Antarctica)) +
  geom_hline(yintercept = rmse_base["East_Antarctica"], col = "red", lty = "dashed")
ggplot(rmse_data) +
  geom_line(aes(x = iter, y = Weddell)) +
  geom_hline(yintercept = rmse_base["Weddell"], col = "red", lty = "dashed")
ggplot(rmse_data) +
  geom_line(aes(x = iter, y = Bellingshausen_Amundsen_Sea)) +
  geom_hline(yintercept = rmse_base["Bellingshausen_Amundsen_Sea"], col = "red", lty = "dashed")


# check results graphically
plot(results$a[[pl[1]]][1,1:(3*tau)], type = "l")
lines(results$a[[pl[2]]][1,1:(3*tau)], col = 2)
lines(results$a[[pl[3]]][1,1:(3*tau)], col = 3)

plot(results$a[[pl[1]]][2,1:(3*tau)], type = "l")
lines(results$a[[pl[2]]][2,1:(3*tau)], col = 2)
lines(results$a[[pl[3]]][2,1:(3*tau)], col = 3)



plot(results$b[[pl[1]]][1,1:(3*tau)], type = "l")
lines(results$b[[pl[2]]][1,1:(3*tau)], col = 2)
lines(results$b[[pl[3]]][1,1:(3*tau)], col = 3)

plot(results$b[[pl[1]]][2,1:(3*tau)], type = "l")
lines(results$b[[pl[2]]][2,1:(3*tau)], col = 2)
lines(results$b[[pl[3]]][2,1:(3*tau)], col = 3)


plot(results$vecH[[pl[1]]], type = "l")
lines(results$vecH[[pl[2]]], col = 2)
lines(results$vecH[[pl[3]]], col = 3)


# compare these two, this seems reasonable
lm_H = lm(a[1,] ~ -1 + t(b))
summary(lm_H)
lm_H
H[1,]



# par(mfrow = c(1, 3))
# image(results$H[[pl[1]]])
# image(results$H[[pl[2]]])
# image(results$H[[pl[3]]])
# par(mfrow = c(1, 1))


plot(results$vecB[[pl[1]]], type = "l")
lines(results$vecB[[pl[2]]], col = 2)
lines(results$vecB[[pl[3]]], col = 3)

par(mfrow = c(1, 3))
image(results$B[[pl[1]]])
image(results$B[[pl[2]]])
image(results$B[[pl[3]]])
par(mfrow = c(1, 1))


# r seems really low?
plot(results$r[-(1:10)], type = "l")
# and s maybe a little high?
plot(results$s[-(1:10)], type = "l")

plot(as.vector(results$C[[pl[1]]]), type = "l")
lines(as.vector(results$C[[pl[2]]]), col = 2)
lines(as.vector(results$C[[pl[3]]]), col = 3)

par(mfrow = c(1, 3))
image(results$C[[pl[1]]])
image(results$C[[pl[2]]])
image(results$C[[pl[3]]])
par(mfrow = c(1, 1))


plot(as.vector(results$D[[pl[1]]]), type = "l")
lines(as.vector(results$D[[pl[2]]]), col = 2)
lines(as.vector(results$D[[pl[3]]]), col = 3)

par(mfrow = c(1, 3))
image(results$D[[pl[1]]])
image(results$D[[pl[2]]])
image(results$D[[pl[3]]])
par(mfrow = c(1, 1))

