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
daily_cov_daily = cov_s

