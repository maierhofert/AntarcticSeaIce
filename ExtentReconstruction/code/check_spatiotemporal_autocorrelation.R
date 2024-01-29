# multivariate version of bk.1.s
library(bayesbackcast)
library(rstan)
library("dplyr")
fname <- c("01", "02", "03", "04", "05")
#  bump up to 10, was 4 and 4 originally
n_cores <- 10
chains <- 10
sel_region <- match(fname,c("01","02","03","04","05"))
fname <- "01.02..all"
vname <- "s"
## To avoid recompilation of unchanged Stan programs, we recommend calling
rstan_options(auto_write = TRUE, cores = n_cores)
options(mc.cores=n_cores)
library(bayesplot)
library(loo)
set.seed(1)

# read in data
# all_dat = readRDS("../data/complete_dataset/all_data.RDS")
# # when running from within the bayesbackcast package
# all_dat = readRDS("../RECONSTRUCTION/data/complete_dataset/all_data.RDS")
all_dat = readRDS("../data/complete_dataset/all_dataLag.RDS")


# subset to before 2021
all_dat = all_dat[all_dat$Year <= 2020,]
# # subset to first 100 PC
# all_dat = all_dat[,1:109]
all_dat$tdate = all_dat$Year + (all_dat$Month -0.5) / 12

# extASN <- grep("ASN0",colnames(all_dat),fixed=T)
# all_dat <- all_dat[,-extASN]
#inc <- colnames(all_dat)[c(1:8,grep("ppb",colnames(all_dat),fixed=TRUE))]
#inc <- c(inc, "VMarsh_O_Higgins_SLP","VMarsh_O_Higgins_tmp","VScott_Base_mcMurdo_SLP","VMarsh_O_Higgins_tmp")
#inc <- inc[-grep("_l7",inc,fixed=T)]
#inc
#all_dat <- all_dat[,inc]
all_x <- all_dat[,-c(1:9)]
#round(cor(all_dat[,inc],use="complete"),2)
#buff.start <- which.max(all_dat$tdate > 1978.82)-1
#buff <- all_dat$tdate >= all_dat$tdate[buff.start] & all_dat$tdate < 2021.5
buff <- all_dat$Date >= as.Date("1979-01-01") & all_dat$Date < as.Date("2021-07-30")
# add a couple more year before
# will we get predictions for it?
buff <- all_dat$Date >= as.Date("1899-01-01") & all_dat$Date < as.Date("2021-07-30")

#odd <- rep(c(TRUE, FALSE), length.out=length(buff))
#buff <- all_dat$tdate >= all_dat$tdate[buff.start] & all_dat$tdate < 2021.5 & odd
all_dat <- all_dat[buff,]
all_x <- all_x[buff,]
all_dat$sie <- all_dat[,sel_region+3]
# how many months are missing
sum(is.na(all_dat$sie)) / 5
# how many years are mssing
sum(is.na(all_dat$sie)) / 5 / 12

# standardize sea ice extent data using monthly sd
sie_dat = all_dat$sie
sie_dat$Month = all_dat$Month
# monthly_mean = sie_dat %>%
#   group_by(Month) %>%
#   summarize_all(.funs = list(mean = mean),na.rm = TRUE)
# monthly_mean_df = do.call("rbind", replicate(122, monthly_mean, simplify = FALSE))
monthly_sd = sie_dat %>%
  group_by(Month) %>%
  summarize_all(.funs = list(sd = sd),na.rm = TRUE)
# plot the standard deviations
matplot(monthly_sd[,-1], type = "l")
monthly_sd_df = do.call("rbind", replicate(122, monthly_sd, simplify = FALSE))
# sea ice is already demeaned
# sie_dat_demean = sie_dat[,1:5] - monthly_mean_df[2:6]
sie_dat_stand = sie_dat[,1:5] / monthly_sd_df[2:6]
summary(sie_dat_stand)
apply(sie_dat_stand, 2, sd, na.rm = TRUE)
all_dat$sie = sie_dat_stand

all_dat$sie[is.na(all_dat$sie)] <- -1000
sie <- ts(data = all_dat$sie, start = min(all_dat$tdate), end = max(all_dat$tdate),
          frequency = 12)

# ##################
all_dat$observed = FALSE
all_dat$observed[958:1464] = TRUE

# ##########################
# check covariances over the course of the year
cov(all_dat[(all_dat$Month == 1) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 2) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 3) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 4) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 5) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 6) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 7) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 8) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 9) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 10) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 11) & all_dat$observed, "sie"]) %>%
  round(2)
cov(all_dat[(all_dat$Month == 12) & all_dat$observed, "sie"]) %>%
  round(2)

spat_covs = list()
for (i in 1:12) {
  spat_covs[[i]] = cov(all_dat[(all_dat$Month == i) & all_dat$observed, "sie"])
}

spat_covs = array(NA, dim = c(5, 5, 12))
for (i in 1:12) {
  spat_covs[,,i] = cov(all_dat[(all_dat$Month == i) & all_dat$observed, "sie"])
}
sd(spat_covs[1,2,])
sd(spat_covs[1,3,])
sd(spat_covs[1,4,])
sd(spat_covs[1,5,])

sd(spat_covs[2,3,])
sd(spat_covs[2,4,])
sd(spat_covs[2,5,])

mean(spat_covs[3,4,])
sd(spat_covs[3,4,])
sd(spat_covs[3,5,])

sd(spat_covs[4,5,])

# ##########################
# check temporal autocorrelation over the course of the year
cov_dat = sie_dat_stand[958:1464,]
cov_dat$Month = c(10:12, rep(1:12, 42))

# January King_Hakon
cor(cov_dat[cov_dat$Month == 1, 1],
    cov_dat[cov_dat$Month == 12, 1][-43])
# February King_Hakon
cor(cov_dat[cov_dat$Month == 2, 1],
    cov_dat[cov_dat$Month == 1, 1])
# March King_Hakon
cor(cov_dat[cov_dat$Month == 3, 1],
    cov_dat[cov_dat$Month == 2, 1])

# correlate with lag 1
cov_dat_l1 = dplyr::lag(cov_dat, 1)
cor(cov_dat[cov_dat$Month == 1, 1],
    cov_dat_l1[cov_dat$Month == 1, 1])

# King Haakon
KH_lag1 = rep(NA, 12)
for (i in 1:12) {
  KH_lag1[i] = cor(cov_dat[cov_dat$Month == i, 1],
                   cov_dat_l1[cov_dat$Month == i, 1], use = "pairwise.complete.obs")
}
mean(KH_lag1)
sd(KH_lag1)

# Ross
Ross_lag1 = rep(NA, 12)
for (i in 1:12) {
  Ross_lag1[i] = cor(cov_dat[cov_dat$Month == i, 2],
                   cov_dat_l1[cov_dat$Month == i, 2], use = "pairwise.complete.obs")
}
mean(Ross_lag1)
sd(Ross_lag1)

# East Antarctica
EA_lag1 = rep(NA, 12)
for (i in 1:12) {
  EA_lag1[i] = cor(cov_dat[cov_dat$Month == i, 3],
                     cov_dat_l1[cov_dat$Month == i, 3], use = "pairwise.complete.obs")
}
mean(EA_lag1)
sd(EA_lag1)

# Weddell
Weddell_lag1 = rep(NA, 12)
for (i in 1:12) {
  Weddell_lag1[i] = cor(cov_dat[cov_dat$Month == i, 4],
                     cov_dat_l1[cov_dat$Month == i, 4], use = "pairwise.complete.obs")
}
mean(Weddell_lag1)
sd(Weddell_lag1)

# BA
BA_lag1 = rep(NA, 12)
for (i in 1:12) {
  BA_lag1[i] = cor(cov_dat[cov_dat$Month == i, 5],
                     cov_dat_l1[cov_dat$Month == i, 5], use = "pairwise.complete.obs")
}
mean(BA_lag1)
sd(BA_lag1)

# ###############################################################################

# correlate with lag 2
cov_dat_l2 = dplyr::lag(cov_dat, 2)
KH_lag2 = rep(NA, 12)
for (i in 1:12) {
  KH_lag2[i] = cor(cov_dat[cov_dat$Month == i, 1],
                   cov_dat_l2[cov_dat$Month == i, 1], use = "pairwise.complete.obs")
}
mean(KH_lag2)
sd(KH_lag2)

# Ross
Ross_lag2 = rep(NA, 12)
for (i in 1:12) {
  Ross_lag2[i] = cor(cov_dat[cov_dat$Month == i, 2],
                     cov_dat_l2[cov_dat$Month == i, 2], use = "pairwise.complete.obs")
}
mean(Ross_lag2)
sd(Ross_lag2)


# East Antarctica
EA_lag2 = rep(NA, 12)
for (i in 1:12) {
  EA_lag2[i] = cor(cov_dat[cov_dat$Month == i, 3],
                   cov_dat_l2[cov_dat$Month == i, 3], use = "pairwise.complete.obs")
}
mean(EA_lag2)
sd(EA_lag2)

# Weddell
Weddell_lag2 = rep(NA, 12)
for (i in 1:12) {
  Weddell_lag2[i] = cor(cov_dat[cov_dat$Month == i, 4],
                        cov_dat_l2[cov_dat$Month == i, 4], use = "pairwise.complete.obs")
}
mean(Weddell_lag2)
sd(Weddell_lag2)

# BA
BA_lag2 = rep(NA, 12)
for (i in 1:12) {
  BA_lag2[i] = cor(cov_dat[cov_dat$Month == i, 5],
                   cov_dat_l2[cov_dat$Month == i, 5], use = "pairwise.complete.obs")
}
mean(BA_lag2)
sd(BA_lag2)



