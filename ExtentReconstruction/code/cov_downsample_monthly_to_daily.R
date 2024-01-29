# read in covariates
amelia_monthly_group6 = readRDS("../data/fogt_predictor_data_imputed/amelia_monthly_group6.RDS")
amelia_monthly_group83 = readRDS("../data/fogt_predictor_data_imputed/amelia_monthly_group83.RDS")
amelia_monthly_group88 = readRDS("../data/fogt_predictor_data_imputed/amelia_monthly_group88.RDS")
amelia_monthly_group9 = readRDS("../data/fogt_predictor_data_imputed/amelia_monthly_group9.RDS")

cov1 = amelia_monthly_group6$imputations[["imp1"]][,c(21, 1:20)]
cov2 = amelia_monthly_group83$imputations[["imp1"]][,3:50]
cov3 = amelia_monthly_group88$imputations[["imp1"]][,3:41]
cov4 = amelia_monthly_group9$imputations[["imp1"]][,3:62]
cov = cbind(cov1, cov2, cov3, cov4)
summary(cov)

# do an eigenvector decomposition on covariates
# maybe do a smoothed version of this
eigen_X = prcomp(cov[-(1:6),-(1:3)])

# use the first ??? principal components, see
plot(eigen_X$sdev)
abline(h = 1)
sum(eigen_X$sdev > 1)

cor_X = cor(cov[-(1:6),-(1:3)])
image(cor_X)
image(abs(cor_X))

# there is still structure in the first 15 or so I think
image(abs(eigen_X$rotation[,1:20]), xlab = "Covariates", ylab = "Eigenvectors")
image(abs(eigen_X$rotation[,1:10]), xlab = "Covariates", ylab = "Eigenvectors")
image(abs(eigen_X$rotation[,11:20]), xlab = "Covariates", ylab = "Eigenvectors")
image(abs(eigen_X$rotation[,1:50]), xlab = "Covariates", ylab = "Eigenvectors")


# smooth interpolation
Dates <- seq.Date(lubridate::ymd("1900-01-01"), lubridate::ymd("2018-12-31"), by = 1)
cov_daily = data.frame(Date = Dates)

# for (i in 4:ncol(cov)) {
#   ###cubic spline for comparison
#   SplineFun <- splinefun(x = cov$Date, y = cov[,i])
#   SplineFit <- SplineFun(Dates)
#   cov_daily[, names(cov)[i]] = SplineFit
# }
for (i in 4:ncol(cov)) {
  smooth = ksmooth(cov$Date, cov[,i], kernel = "normal", bandwidth = 30.4)$y
  SplineFun <- splinefun(x = cov$Date, y = smooth)
  SplineFit <- SplineFun(Dates)
  cov_daily[, names(cov)[i]] = SplineFit
}
head(cov_daily)

plot(x = cov$Date, y = cov[,4], xlim = c(lubridate::ymd("1901-01-01"), lubridate::ymd("1910-01-01")))
lines(cov_daily$Date, cov_daily[,4 - 2], col = "red")
# lines(cov$Date, smooth, col = "blue")
apply(cov_daily, 2, sd)

mon_cov_daily = cov_daily
# # ###############
# # the PC does not change by interpolation
# # do an eigenvector decomposition on covariates
# # maybe do a smoothed version of this
# eigen_X = prcomp(cov_daily[-(1:20000),-(1:3)])
# 
# # use the first ??? principal components, see
# plot(eigen_X$sdev)
# abline(h = 1)
# sum(eigen_X$sdev > 1)
# 
# cor_X = cor(cov_daily[-(1:20000),-(1:3)])
# image(cor_X)
# image(abs(cor_X))




