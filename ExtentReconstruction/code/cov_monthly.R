# read in covariates
amelia_monthly_group6 = readRDS("../data/fogt_predictor_data_imputed/amelia_monthly_group6.RDS")
amelia_monthly_group8 = readRDS("../data/fogt_predictor_data_imputed/amelia_monthly_group8.RDS")
amelia_monthly_group9 = readRDS("../data/fogt_predictor_data_imputed/amelia_monthly_group9.RDS")

cov1 = amelia_monthly_group6$imputations[["imp1"]][,c(21, 1:20)]
cov2 = amelia_monthly_group8$imputations[["imp1"]][,3:84]
cov3 = amelia_monthly_group9$imputations[["imp1"]][,3:66]
cov = cbind(cov1, cov2, cov3)
cov = cov[,c(1:3, order(names(cov)[-c(1:3)]) + 3)]
summary(cov)
# standardize variables
cov = cov[order(cov$Date),]
# check out sd
sds = apply(cov[,-(1:3)], 2, sd)
hist(sds)
# check out means
means = colMeans(cov[,-(1:3)])
hist(means)

# standardize
cov2 = cov
for (i in 4:ncol(cov)) {
  cov2[,i] = (cov[,i] - means[names(cov)[i]]) / sds[names(cov)[i]]
}
cov = cov2
# save covariates
saveRDS(cov, "../data/fogt_predictor_data_imputed/all_monthly_cov.RDS")


# do an eigenvector decomposition on covariates
# maybe do a smoothed version of this
eigen_cov = prcomp(cov[,-(1:3)])
# use the first ??? principal components, see
plot(eigen_cov$sdev)
abline(h = 1)
sum(eigen_cov$sdev > 1)

# compute correlation matrix
cor_cov = cor(cov[,-(1:3)])
# plot correlation matrix
library(ggplot2)
melted_cor_cov <- reshape2::melt(cor_cov)
ggplot(data = melted_cor_cov, aes(x=Var1, y=Var2, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Correlation", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../plots/covariates/all_covariates/corr.png", device = "png", width = 25, height = 25)
ggplot(data = melted_cor_cov, aes(x=Var1, y=Var2, fill=abs(value))) + 
  xlab("") + ylab("") +
  scale_fill_gradient("abs(Correlation)", limits = c(0, 1), 
                      high = "firebrick", low = "white") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../plots/covariates/all_covariates/abs_corr.png", device = "png", width = 25, height = 25)

# plot the eigenvectors
melted_eigenvec <- reshape2::melt(eigen_cov$rotation)
ggplot(data = melted_eigenvec, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() +
  xlim(paste0("PC", 1:75)) +
  ylab("") + ylab("") +
  scale_fill_gradient2("Loading", high = "firebrick", low = "steelblue", limits = c(-0.25, 0.25), oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../plots/covariates/all_covariates/eigenvec_col.png", device = "png", width = 25, height = 25)

ggplot(data = melted_eigenvec, aes(x=Var2, y=Var1, fill=abs(value))) + 
  geom_tile() +
  xlim(paste0("PC", 1:75)) +
  ylab("") + ylab("") +
  scale_fill_gradient("abs(Loading)", high = "black", low = "white", limits = c(0, 0.25), oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../plots/covariates/all_covariates/eigenvec_bw.png", device = "png", width = 25, height = 25)

# ###############################################################################
# for (i in 4:ncol(cov)) {
#   # for (i in 4:10) {
#   pacf(cov[,i], main = names(cov)[i], lag.max = 13)
# }


# add lagged versions of all covariates
cov = readRDS("../data/fogt_predictor_data_imputed/all_monthly_cov.RDS")

lags = c(1, 2, 3, 4, 12)
col_names = names(cov)
jmax = ncol(cov)
# save_date = cov$Date
# add buffer of random variables for 1898
cov = rbind(matrix(rnorm(12 * jmax), 
                   ncol = jmax, nrow = 12,
                   dimnames = list(NULL, colnames(cov))),
            cov)
for (j in 4:jmax) {
  for (i in seq_along(lags)) {
    cov[,paste0(col_names[j], "_l", lags[i])] = lag(cov[,col_names[j]], n = lags[i])
  }
}
cov = cov[-(1:12),]
cov$Date = as.Date(paste0(cov$Year, "-", cov$Month, "-15"))
saveRDS(cov, "../data/fogt_predictor_data_imputed/all_monthly_cov_lagged.RDS")

# cov$Date = save_date
# compute PCA
eigen_cov = prcomp(cov[,-(1:3)])
eigen_x = predict(eigen_cov, newdata = cov)
covPC = data.frame(Date = cov$Date, Year = cov$Year, Month = cov$Month, eigen_x)
saveRDS(covPC, "../data/fogt_predictor_data_imputed/all_monthly_PC_lagged.RDS")


melted_eigenvec <- reshape2::melt(eigen_cov$rotation)

# plot eigenvectors
ggplot(data = melted_eigenvec, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() +
  # xlim(paste0("PC", 1:75)) +
  # xlim(paste0("PC", 1:50)) +
  # ylim(names(cov)[50:90]) +
  ylab("") + xlab("") +
  scale_fill_gradient2("Loading", high = "firebrick", low = "steelblue", limits = c(-0.25, 0.25), oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../plots/covariates/all_covariates/eigenvec_lags_col.png", device = "png", width = 25, height = 100, limitsize = FALSE)



