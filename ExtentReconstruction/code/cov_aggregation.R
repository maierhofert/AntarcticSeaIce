# analyze the combined covariates available
source("cov_daily.R")
source("cov_downsample_monthly_to_daily.R")
# TODO there is a SST reconstruction that I am not using currently
# source("cov_climate_indices.R")

# combine all covariates
cov = full_join(daily_cov_daily, mon_cov_daily, by = "Date")
# cov = full_join(cov, ind_cov_daily, by = "Date")
cov = cov[order(cov$Date),]
summary(is.na(cov))

# 2019 is not available for monthly data, but for daily and for indices
# we also do not have the sea ice in 2019, so that should be alright for now
cov$Date[which(is.na(cov$'619010_SLP'))]
cov = cov[cov$Date < as.Date("2019-01-01"),]

# plot all covariates
# some are way wigglier than others
for (i in 2:ncol(cov)) {
  plot(cov$Date, cov[, i], type = "l", main = names(cov)[i])
}

# check out sd
sds = apply(cov[,-1], 2, sd)
hist(sds)
# check out means
means = colMeans(cov[,-1])
hist(means)

# standardize
cov2 = cov
for (i in 2:ncol(cov)) {
  cov2[,i] = (cov[,i] - means[i-1]) / sds[i-1]
}
cov = cov2


# do an eigenvector decomposition on covariates
# maybe do a smoothed version of this
eigen_cov = prcomp(cov[,-1])
# use the first ??? principal components, see
plot(eigen_cov$sdev)
abline(h = 1)
sum(eigen_cov$sdev > 1)

# compute correlation matrix
cor_cov = cor(cov[,-1])
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
# # take out large scale climate indices
# cov = cov[,-c(256, 270)]
# eigen_cov = prcomp(cov[,-1])

# add lagged versions of all covariates
lags = c(2, 7, 29, 61)
col_names = names(cov)
jmax = ncol(cov)
save_date = cov$Date
cov = rbind(matrix(rnorm(61 * jmax), 
                   ncol = jmax, nrow = 61,
                   dimnames = list(NULL, colnames(cov))),
            cov)
for (j in 2:jmax) {
  for (i in seq_along(lags)) {
    cov[,paste0(col_names[j], "_l", lags[i])] = lag(cov[,col_names[j]], n = lags[i])
  }
}
cov = cov[-(1:61),]
cov$Date = save_date
# compute PCA on every 10th observation to speed up computation
eigen_cov = prcomp((cov[,-1][c(T, rep(F, 9)),]))
eigen_x = predict(eigen_cov, newdata = cov)

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

