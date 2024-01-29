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


# get data
#start.date <- which.max(all_dat$tdate > 1978.82)
#start.date
# start <- all_dat$tdate > 1979.05
# all_dat$tdate[start][1:10]
# no1978 <- seq_along(sie[,1])[start]
# table(all_dat$year[no1978])
# all_x = all_x[no1978,]
# sie[all_dat$tdate > 1987.91 & all_dat$tdate < 1988] <- -1000
# sie[all_dat$tdate > 1987.55 & all_dat$tdate < 1987.7] <- -1000
# sum(sie[no1978] < -999)


# # TM: Add splines
# p <- ncol(all_x) - 1
# p
# call_x = colnames(all_x)[1:p]
# all_x_spline = cbind(all_x[1:p],all_x[1:p],all_x[1:p],all_x[1:p])
# Month <- rep(1:12, nrow(all_x)/12)
# pmsMat = mgcv::cSplineDes(x = 1:12, knots = c(0.5, 4.5, 8.5, 12.5), ord = 3, derivs=0)
# pmsMat = mgcv::cSplineDes(x = 1:12, knots = c(0.5, 3.5, 6.5, 9.5, 12.5), ord = 4, derivs=0)
# matplot(pmsMat, type = "l")
# for(j in 1:4){
#   colnames(all_x_spline)[(1:p)+(j-1)*p] <- paste0("b",j,"_",call_x)
#   all_x_spline[,(1:p)+(j-1)*p] <- sweep(all_x_spline[,(1:p)+(j-1)*p],1,pmsMat[Month,j],"*")
# }
# all_x = all_x_spline


# the random forest as a comparisomn
# # A tibble: 5 × 3
# sector                      rmse0  rmse
# <fct>                       <dbl> <dbl>
# 2 Ross                        0.346 0.317
# 5 Bellingshausen_Amundsen_Sea 0.171 0.154
# 4 Weddell                     0.307 0.256
# 1 King_Hakon                  0.294 0.277
# 3 East_Antarctica             0.178 0.169


p <- ncol(all_x)
n <- nrow(all_x)

# These horseshoe parameters without seasonal effect
p_nonzero <- 15 # prior guess for the number of relevant variables
slab_scale = 0.025
# result in res_table
# >   res_table[c(2, 5, 4, 1, 3, 6),]
#               Sector cor_clim cor_mod ce_clim ce_mod rmse_clim rmse_mod
# 2           Ross Sea        0    0.34       0   0.10      0.34     0.33
# 5 Amundsen Bellings.        0    0.53       0   0.27      0.17     0.14
# 4        Weddell Sea        0    0.61       0   0.37      0.30     0.24
# 1    King Haakon VII        0    0.33       0   0.10      0.29     0.28
# 3    East Antarctica        0    0.31       0   0.06      0.18     0.17
# 6              Total        0    0.09       0  -0.06      0.56     0.58
# >   this.loo
#
# Computed from 2500 by 504 log-likelihood matrix
#
# Estimate   SE
# elpd_loo   1314.1 48.5
# p_loo       277.0 10.5
# looic     -2628.3 97.1
# ------
#   Monte Carlo SE of elpd_loo is NA.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     456   90.5%   5
# (0.5, 0.7]   (ok)        41    8.1%   3
# (0.7, 1]   (bad)        7    1.4%   1
# (1, Inf)   (very bad)   0    0.0%   <NA>
#   See help('pareto-k-diagnostic') for details.

# These horseshoe parameters without seasonal effect
# This is currently in the paper
p_nonzero <- 15
slab_scale = 0.25
# >   res_table[c(2, 5, 4, 1, 3, 6),]
# Sector cor_clim cor_mod ce_clim ce_mod rmse_clim rmse_mod
# 2           Ross Sea        0    0.28       0   0.04      0.34     0.34
# 5 Amundsen Bellings.        0    0.54       0   0.28      0.17     0.14
# 4        Weddell Sea        0    0.62       0   0.38      0.30     0.24
# 1    King Haakon VII        0    0.38       0   0.14      0.29     0.27
# 3    East Antarctica        0    0.32       0   0.07      0.18     0.17
# 6              Total        0    0.05       0  -0.11      0.56     0.59

# >   this.loo
#
# Computed from 2500 by 504 log-likelihood matrix
#
# Estimate   SE
# elpd_loo   1325.1 48.4
# p_loo       282.4 10.6
# looic     -2650.2 96.9
# ------
#   Monte Carlo SE of elpd_loo is NA.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     448   88.9%   5
# (0.5, 0.7]   (ok)        51   10.1%   2
# (0.7, 1]   (bad)        5    1.0%   2
# (1, Inf)   (very bad)   0    0.0%   <NA>
#   See help('pareto-k-diagnostic') for details.

# # if we leave out the MA component, strange inclonclusive differences
# Sector cor_clim cor_mod ce_clim ce_mod rmse_clim rmse_mod
# 2           Ross Sea        0    0.23       0   0.00      0.34     0.34
# 5 Amundsen Bellings.        0    0.58       0   0.33      0.17     0.14
# 4        Weddell Sea        0    0.63       0   0.39      0.30     0.24
# 1    King Haakon VII        0    0.36       0   0.13      0.29     0.27
# 3    East Antarctica        0    0.34       0   0.09      0.18     0.17
# 6              Total        0    0.05       0  -0.12      0.56     0.59

# These horseshoe parameters with seasonal effect
p_nonzero <- 15*4
slab_scale = 0.025
# >   res_table[c(2, 5, 4, 1, 3, 6),]
# Sector cor_clim cor_mod ce_clim ce_mod rmse_clim rmse_mod
# 2           Ross Sea        0    0.18       0   0.03      0.34     0.34
# 5 Amundsen Bellings.        0    0.33       0   0.08      0.17     0.16
# 4        Weddell Sea        0    0.40       0   0.15      0.30     0.28
# 1    King Haakon VII        0    0.12       0   0.01      0.29     0.29
# 3    East Antarctica        0    0.32       0   0.09      0.18     0.17
# 6              Total        0    0.04       0  -0.06      0.56     0.58
# >   this.loo
#
# Computed from 2500 by 504 log-likelihood matrix
#
# Estimate   SE
# elpd_loo   1232.4 45.5
# p_loo       196.4 10.0
# looic     -2464.7 91.0
# ------
#   Monte Carlo SE of elpd_loo is NA.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     436   86.5%   5
# (0.5, 0.7]   (ok)        60   11.9%   2
# (0.7, 1]   (bad)        7    1.4%   1
# (1, Inf)   (very bad)   1    0.2%   1
# See help('pareto-k-diagnostic') for details.

p_nonzero <- 30*4
lab_scale = 0.025
# >   res_table[c(2, 5, 4, 1, 3, 6),]
# Sector cor_clim cor_mod ce_clim ce_mod rmse_clim rmse_mod
# 2           Ross Sea        0    0.20       0   0.04      0.34     0.34
# 5 Amundsen Bellings.        0    0.36       0   0.11      0.17     0.16
# 4        Weddell Sea        0    0.40       0   0.15      0.30     0.28
# 1    King Haakon VII        0    0.11       0   0.00      0.29     0.29
# 3    East Antarctica        0    0.25       0   0.03      0.18     0.17
# 6              Total        0    0.06       0  -0.05      0.56     0.57
# >   this.loo
#
# Computed from 2500 by 504 log-likelihood matrix
#
# Estimate   SE
# elpd_loo   1229.2 45.8
# p_loo       199.0 10.6
# looic     -2458.4 91.6
# ------
#   Monte Carlo SE of elpd_loo is NA.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     435   86.3%   6
# (0.5, 0.7]   (ok)        57   11.3%   2
# (0.7, 1]   (bad)       10    2.0%   1
# (1, Inf)   (very bad)   2    0.4%   0
# See help('pareto-k-diagnostic') for details.

# number of autoregressive components
Kar <- 4
# 0 moving average components
Kma <- 0 # no MA component
mdl <- arima(ts = sie, order = c(Kar,0,Kma),
             xreg = as.matrix(all_x[,1:984]), # all_x[,1:984], # xreg = all_pc[,1:50],
             slab_scale=slab_scale,
             tprior=FALSE,
             series.name = "all", seldat = TRUE, #seldat=no1978,
             lagAR=c(1,2,3,12),lagMA=0, p_nonzero=p_nonzero)
str(mdl)

# # this recompiles stan code
# pkgbuild::clean_dll(path = ".")
# pkgbuild::compile_dll(path = ".", force = T)
# then hit build -> more -> load all
# then hit build -> install and restart

# # not sure about any of this
# # should check and only recompile new files
# pkgbuild::compile_dll(path = ".", force = F)
# # or
# pkgbuild::compile_dll(path = ".", force = F)
# devtools::install(quick=FALSE)

# mdl = mdl[names(mdl) != "n"]
# attr(mdl, "class") = "arima"
# str(mdl)
# mdl <- reg(sie,xreg = all_pc #slab_scale=0.025,
#           series.name = "all",seldat=no1978, p_nonzero=p_nonzero)

# # run some tests
# mdl$yall = mdl$yall[-(1:6),]
# mdl$N = 10L
# mdl$test=1
if(T){
  # varstan calls  # bayesbackcast:::fit_arima
  # bayesbackcast:::stanmodels$arima_nostd
  # this has a 26 minute runtime for all cov with lags
  fit <- varstan(model = mdl, iter=20, chains=1, tree.depth=10)
  # this has a 5h runtime
  # before iter=1000
  fit <- varstan(mode = mdl, iter=500, chains=chains, tree.depth=10)
  # saveRDS(fit, file = paste0("../output/bc.",fname,".od.",vname,".rds"))
  # # save location when in bayesbackcast package
  # saveRDS(fit, file = paste0("../RECONSTRUCTION/output/bc.",fname,".od.",vname, ".nz.", p_nonzero, "IncMissingAllCovLag.rds"))
  # # save location when in RECONSTRUCTION/code directory
  # saveRDS(fit, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLag.rds"))
  # saveRDS(fit, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagSeasonal.rds"))
  # saveRDS(fit, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagSeasonalSeasonStand.rds"))
  saveRDS(fit, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagSeasonalSeasonStandSeasonalAR.rds"))
  # saveRDS(fit, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagSeasonStand.rds"))
} else{
  fit <- readRDS(file = paste0("../RECONSTRUCTION/output/bc.",fname,".od.",vname, ".nz.", p_nonzero, "IncMissingAllCovLag.rds"))
  # very little modeled variation
  fit = readRDS(file = paste0("../output/bc.01.02..all.od.s.nz.60.ss.0.025IncMissingAllCovLagSeasonalSeasonStand.rds"))
  # this is currently in paper, no seasonally varying coefs
  fit = readRDS(file = paste0("../output/bc.01.02..all.od.s.nz.15.ss.0.25IncMissingAllCovLagSeasonStand.rds"))
  # seasonally varying ar coefficients
  fit = readRDS(file = paste0("../output/bc.01.02..all.od.s.nz.15.ss.0.25IncMissingAllCovLagSeasonalSeasonStandSeasonalAR.rds"))
}

# do a cross validation
library("foreach")
if (F) {
  # register a parallel backend
  ncores = 20
  ncores = 10

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  yrs = 1979:2020

  # cross validate climatology as well
  clim_avg_cv = all_dat[961:1464,c(4:8, 3)]
  sie_obs = clim_avg_cv
  for (i in 1:nrow(sie_obs)) {
    month_ind = (1:nrow(sie_obs)) %% 12
    this.month_ind = i %% 12
    clim_avg_cv[i,] = colMeans(sie_obs[(month_ind == this.month_ind) & !((1:nrow(sie_obs)) %in% c(i-2, i-1, i, i+1, i+2)),])
  }
  # colMeans(abs(clim_avg_cv))



  CVfits = NA
  CVfits = foreach (i = 1:length(yrs), .packages = "bayesbackcast") %dopar% {
    # CVfits = foreach (i = 1:10, .packages = "bayesbackcast") %dopar% {
    del_yrs = (-2):2 + yrs[i]
    this.sie = sie
    this.sie[all_dat$Year %in% del_yrs,] = -1000
    this.sie = window(this.sie, start = 1979)
    this.all_x = all_x[all_dat$Year >= 1979,]
    n = nrow(this.sie)
    this.mdl <- bayesbackcast::arima(ts = this.sie,
                                     order = c(Kar,0,Kma),
                                     xreg = as.matrix(this.all_x[,1:984]), # this.all_pc[,1:50],
                                     slab_scale=slab_scale,
                                     tprior=FALSE,
                                     series.name = "all",
                                     seldat = TRUE, #seldat=no1978,
                                     lagAR=c(1,2,3,12),lagMA=1, p_nonzero=p_nonzero)
    # check out traceplots
    this.fit <- varstan(mode = this.mdl, iter=400, chains=1, tree.depth=10)
    return(this.fit)
  }
  # saveRDS(CVfits, file = paste0("../RECONSTRUCTION/output/bc.",fname,".od.",vname, ".nz.", p_nonzero, "IncMissingAllCovLagCVfitsSeasonal.rds"))
  # saveRDS(CVfits, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagCVfitsSeasonal.rds"))
  # saveRDS(CVfits, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagCVfitsSeasonStand.rds"))
  saveRDS(CVfits, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagCVfitsSeasonalSeasonStandSeasonalAr.rds"))
  # saveRDS(CVfits, file = paste0("../output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".ss.", slab_scale, "IncMissingAllCovLagCVfits.rds"))
  # CVfits = readRDS(file = paste0("../output/bc.01.02..all.od.s.nz.15IncMissingAllCovLagCVfits.rds"))
  # CVfits = readRDS(file = paste0("../output/bc.01.02..all.od.s.nz.60.ss.0.025IncMissingAllCovLagCVfitsSeasonalSeasonStand.rds"))
  #CVfits = readRDS(file = paste0("../output/bc.01.02..all.od.s.nz.15.ss.0.25IncMissingAllCovLagCVfitsSeasonStand.rds"))
  CVfits = readRDS(file = paste0("../output/bc.01.02..all.od.s.nz.15.ss.0.25IncMissingAllCovLagCVfitsSeasonalSeasonStandSeasonalAr.rds"))
  # str(CVfits,1)
  # get predicted values
  preds = lapply(CVfits, function(this.fit) {
    rstan::extract(this.fit$stanfit, "fit")[[1]] # Y_latent
  })
  # combine predictions
  only_predicted = rstan::extract(CVfits[[1]]$stanfit, "fit")[[1]]
  only_predicted[,,] = NA
  this.all_dat = all_dat[all_dat$Year >= 1979,]
  yrs = 1979:2020
  for (i in 1:length(CVfits)) {
    only_predicted[,this.all_dat$Year == yrs[i],] = preds[[i]][,this.all_dat$Year == yrs[i],]
  }
  dim(only_predicted)
  # ###############################################################
  # retransform predictions to be on original scale
  # this is only necessary if they were standardized by month
  # multiply by monthly standard deviation
  for (i in 1:(dim(only_predicted)[1])) {
    only_predicted[i,,] = only_predicted[i,,] * as.matrix(monthly_sd_df[961:1464,2:6])
  }
  # ###################################################
  only_predicted[200,,] %>%
    cov() %>%
    cov2cor() %>%
    round(2)
  # extract(CVfits[[1]]$stanfit, "fit")[[1]][200,,] %>%
  #   cov() %>%
  #   cov2cor() %>%
  #   round(2)
  # extract(CVfits[[1]]$stanfit, "Y_latent")[[1]][200,,] %>%
  #   cov() %>%
  #   cov2cor() %>%
  #   round(2)
  # extract(CVfits[[1]]$stanfit, "epsilon")[[1]][200,,] %>%
  #   cov() %>%
  #   cov2cor() %>%
  #   round(2)

  # # are my predictions different from ground truth?
  # only_predicted[,,1] %>% apply(1,sd) %>% hist()
  # only_predicted[,,1] %>% apply(1,sd) %>% mean()
  # all_dat$King_Hakon %>% sd(na.rm = TRUE)
  #
  # only_predicted[,,2] %>% apply(1,sd) %>% hist()
  # only_predicted[,,2] %>% apply(1,sd) %>% mean()
  # all_dat$Ross %>% sd(na.rm = TRUE)
  #
  # only_predicted[,,3] %>% apply(1,sd) %>% hist()
  # only_predicted[,,3] %>% apply(1,sd) %>% mean()
  # all_dat$East_Antarctica %>% sd(na.rm = TRUE)
  #
  # only_predicted[,,4] %>% apply(1,sd) %>% hist()
  # only_predicted[,,4] %>% apply(1,sd) %>% mean()
  # all_dat$Weddell %>% sd(na.rm = TRUE)
  #
  # only_predicted[,,5] %>% apply(1,sd) %>% hist()
  # only_predicted[,,5] %>% apply(1,sd) %>% mean()
  # all_dat$Bellingshausen_Amundsen_Sea %>% sd(na.rm = TRUE)

  # compute total fit
  total_fits_val = matrix(NA, nrow = dim(only_predicted)[1], ncol = dim(only_predicted)[2])
  for (i in 1:dim(only_predicted)[1]) {
    total_fits_val[i,] = rowSums(only_predicted[i,,])
  }

  # compute average prediction by sector and month
  avg_prediction = data.frame("King_Haakon" = NA, "Ross" = NA,
                              "East_Antarctica" = NA, "Weddell" = NA,
                              "Bellingshausen_Amundsen" = NA)
  for (i in 1:dim(only_predicted)[2]) {
    avg_prediction[i, 1:5] = colMeans(only_predicted[,i,])
  }
  # get a prediction for the total
  avg_prediction$Total = rowSums(avg_prediction)
  # avg_prediction$Total2 = colMeans(total_fits_val)


  # recompute the total as the sum of the sectors
  # the sum of the annual cycles by sector is not the same as the sum of the total annual cycle
  # due to smoothness
  this.all_dat$total = rowSums(this.all_dat[,4:8])

  # plot the out of sample estimates, average, and ground truth
  png(file = "../plots/BayesianModel/LOOCV_predictions%02d.png", width = 7*60, height = 5*60)
  # King Haakon
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,1]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "King Haakon VII", xlab = "Date", ylab = "Anomaly")
  lines(all_dat$tdate[961:1464], avg_prediction[,1], col = 2)
  lines(all_dat$tdate[961:1464], all_dat$King_Hakon[961:1464], col = "steelblue1")
  # Ross
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,2]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Ross Sea", xlab = "Date", ylab = "Anomaly")
  lines(all_dat$tdate[961:1464], avg_prediction[,2], col = 2)
  lines(all_dat$tdate[961:1464], all_dat$Ross[961:1464], col = 4)
  # East Antarctica
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,3]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "East Antarctica", xlab = "Date", ylab = "Anomaly")
  lines(all_dat$tdate[961:1464], avg_prediction[,3], col = 2)
  lines(all_dat$tdate[961:1464], all_dat$East_Antarctica[961:1464], col = "steelblue1")
  # Weddell
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,4]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Weddell Sea", xlab = "Date", ylab = "Anomaly")
  lines(all_dat$tdate[961:1464], avg_prediction[,4], col = 2)
  lines(all_dat$tdate[961:1464], all_dat$Weddell[961:1464], col = "steelblue1")
  # Bellingshausen Amundsen
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,5]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Bellingshausen Amundsen Sea", xlab = "Date", ylab = "Anomaly")
  lines(all_dat$tdate[961:1464], avg_prediction[,5], col = 2)
  lines(all_dat$tdate[961:1464], all_dat$Bellingshausen_Amundsen_Sea[961:1464], col = "steelblue1")
  # Total
  matplot(x = all_dat$tdate[961:1464], y = t(total_fits_val), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Total", xlab = "Date", ylab = "Anomaly")
  lines(all_dat$tdate[961:1464], colMeans(total_fits_val), col = 2)
  lines(all_dat$tdate[961:1464], all_dat$total[961:1464], col = "steelblue1")
  dev.off()
  
  
  # plot the coverage of 95% PI
  png(file = "../plots/BayesianModel/LOOCV_coverage95PI%02d.png", width = 7*60, height = 5*60)
  # coverage of 95% PI
  # King Haakon
  KHPI = t(apply(only_predicted[,,1], 2, quantile, c(0.025, 0.975)))
  inKHPI =  (all_dat$King_Hakon[961:1464] > KHPI[,1]) &
    (all_dat$King_Hakon[961:1464] < KHPI[,2])
  mean(inKHPI)
  # plot
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,1]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "King Haakon VII", xlab = "Date", ylab = "Anomaly")
  matplot(x = all_dat$tdate[961:1464], y = KHPI, type = "l", lty = 2,
          col = "green", add = TRUE)
  lines(all_dat$tdate[961:1464], all_dat$King_Hakon[961:1464], col = "steelblue1")
  points(all_dat$tdate[961:1464][!inKHPI], all_dat$King_Hakon[961:1464][!inKHPI], pch = 16, col = "red")
  
  # Ross
  RossPI = t(apply(only_predicted[,,2], 2, quantile, c(0.025, 0.975)))
  inRossPI =  (all_dat$Ross[961:1464] > RossPI[,1]) &
    (all_dat$Ross[961:1464] < RossPI[,2])
  mean(inRossPI)
  # plot
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,2]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Ross Sea", xlab = "Date", ylab = "Anomaly")
  matplot(x = all_dat$tdate[961:1464], y = RossPI, type = "l", lty = 2,
          col = "green", add = TRUE)
  lines(all_dat$tdate[961:1464], all_dat$Ross[961:1464], col = "steelblue1")
  points(all_dat$tdate[961:1464][!inRossPI], all_dat$Ross[961:1464][!inRossPI], pch = 16, col = "red")
  
  # East Antarctica
  EAPI = t(apply(only_predicted[,,3], 2, quantile, c(0.025, 0.975)))
  inEAPI =  (all_dat$East_Antarctica[961:1464] > EAPI[,1]) &
    (all_dat$East_Antarctica[961:1464] < EAPI[,2])
  mean(inEAPI)
  # plot
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,3]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "East Antarctica", xlab = "Date", ylab = "Anomaly")
  matplot(x = all_dat$tdate[961:1464], y = EAPI, type = "l", lty = 2,
          col = "green", add = TRUE)
  lines(all_dat$tdate[961:1464], all_dat$East_Antarctica[961:1464], col = "steelblue1")
  points(all_dat$tdate[961:1464][!inEAPI], all_dat$East_Antarctica[961:1464][!inEAPI], pch = 16, col = "red")
  
  # Weddell Sea
  WeddellPI = t(apply(only_predicted[,,4], 2, quantile, c(0.025, 0.975)))
  inWeddellPI =  (all_dat$Weddell[961:1464] > WeddellPI[,1]) &
    (all_dat$Weddell[961:1464] < WeddellPI[,2])
  mean(inWeddellPI)
  # plot
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,4]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Weddell Sea", xlab = "Date", ylab = "Anomaly")
  matplot(x = all_dat$tdate[961:1464], y = WeddellPI, type = "l", lty = 2,
          col = "green", add = TRUE)
  lines(all_dat$tdate[961:1464], all_dat$Weddell[961:1464], col = "steelblue1")
  points(all_dat$tdate[961:1464][!inWeddellPI], all_dat$Weddell[961:1464][!inWeddellPI], pch = 16, col = "red")
  
  # Bellingshausen Amundsen
  BAPI = t(apply(only_predicted[,,5], 2, quantile, c(0.025, 0.975)))
  inBAPI =  (all_dat$Bellingshausen_Amundsen_Sea[961:1464] > BAPI[,1]) &
    (all_dat$Bellingshausen_Amundsen_Sea[961:1464] < BAPI[,2])
  mean(inBAPI)
  # plot
  matplot(x = all_dat$tdate[961:1464], y = t(only_predicted[,,5]), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Bellingshausen Amundsen Sea", xlab = "Date", ylab = "Anomaly")
  matplot(x = all_dat$tdate[961:1464], y = BAPI, type = "l", lty = 2,
          col = "green", add = TRUE)
  lines(all_dat$tdate[961:1464], all_dat$Bellingshausen_Amundsen_Sea[961:1464], col = "steelblue1")
  points(all_dat$tdate[961:1464][!inBAPI], all_dat$Bellingshausen_Amundsen_Sea[961:1464][!inBAPI], pch = 16, col = "red")
  
  # Total
  totalPI = t(apply(total_fits_val, 2, quantile, c(0.025, 0.975)))
  intotalPI =  (all_dat$total[961:1464] > totalPI[,1]) &
    (all_dat$total[961:1464] < totalPI[,2])
  mean(intotalPI)
  # plot
  matplot(x = all_dat$tdate[961:1464], y = t(total_fits_val), type = "l", lty = 1,
          col = alpha(1, 0.03), main = "Total", xlab = "Date", ylab = "Anomaly")
  matplot(x = all_dat$tdate[961:1464], y = totalPI, type = "l", lty = 2,
          col = "green", add = TRUE)
  lines(all_dat$tdate[961:1464], all_dat$total[961:1464], col = "steelblue1")
  points(all_dat$tdate[961:1464][!intotalPI], all_dat$total[961:1464][!intotalPI], pch = 16, col = "red")
  
  dev.off()
  

  # # get correlation R and CE per Month and sector
  # # correlation by sector
  # cor(avg_prediction[,1], this.all_dat[,4])
  # cor(avg_prediction[,2], this.all_dat[,5])
  # cor(avg_prediction[,3], this.all_dat[,6])
  # cor(avg_prediction[,4], this.all_dat[,7])
  # cor(avg_prediction[,5], this.all_dat[,8])
  # cor(avg_prediction[,6], this.all_dat[,3])


  # CE by sector
  # helper function, coefficient of efficiency
  my_CE = function(truth, pred) {
    1 - sum((truth - pred)^2) / sum((truth - mean(truth))^2)
  }
  # my_CE(this.all_dat[,4], avg_prediction[,1])
  # my_CE(this.all_dat[,5], avg_prediction[,2])
  # my_CE(this.all_dat[,6], avg_prediction[,3])
  # my_CE(this.all_dat[,7], avg_prediction[,4])
  # my_CE(this.all_dat[,8], avg_prediction[,5])
  # my_CE(this.all_dat[,3], avg_prediction[,6])
  # # my_CE(rowSums(this.all_dat[,4:8]), avg_prediction[,6])

  # get RMSE per Month and sector
  res_data2 = data.frame(
    # model
    "mse_King_Haakon" = rep(NA, dim(only_predicted)[2]),
    "mse_Ross" = NA, "mse_East_Antarctica" = NA,
    "mse_Weddell" = NA, "mse_Bellingshausen_Amundsen_Sea" = NA,
    "mse_Total" = NA,
    # climatology
    "mse0_King_Haakon" = NA, "mse0_Ross" = NA, "mse0_East_Antarctica" = NA,
    "mse0_Weddell" = NA, "mse0_Bellingshausen_Amundsen_Sea" = NA,
    "mse0_Total" = NA)
  for (i in 1:dim(only_predicted)[2]) {
    # res_data2[i, 1:6] = abs(avg_prediction[i,] - this.all_dat[i,c(4:8, 3)])
    # res_data2[i,7:12] = abs(this.all_dat[i,c(4:8, 3)])
    res_data2[i, 1:6] = (avg_prediction[i,] - this.all_dat[i,c(4:8, 3)])^2
    res_data2[i,7:12] = (clim_avg_cv[i,] - this.all_dat[i,c(4:8, 3)])^2 # this.all_dat[i,c(4:8, 3)]^2
  }
  dim(res_data2)
  mse_dat = res_data2 %>% summarize_all(list(mean = mean))
  rmse_dat = sqrt(mse_dat)
  # rmse_dat = sqrt(colMeans(res_data2))

  # create absolute error data
  res_data_abs = sqrt(res_data2)
  ae_dat = res_data_abs %>% summarize_all(list(mean = mean))

  # # covariance of absolute errors
  # cov_ae = res_data_abs[1:6] %>%
  #   cov()
  # colnames(cov_ae) = rownames(cov_ae) = c("King H. VII", "Ross Sea", "East Ant.", "Weddell Sea", "B.A. Sea", "Total")
  # round(cov_ae * 10, 2)
  # print(xtable::xtable(cov_ae[c(2, 5, 4, 1, 3, 6),c(2, 5, 4, 1, 3, 6)] * 10))

  # # covariance of squared errors
  # cov_sqe = res_data2[1:6] %>%
  #   cov()
  # round(cov_sqe * 10, 2)

  # covariance of prediction errors
  res_data_err = data.frame(
    # model
    "mse_King_Haakon" = rep(NA, dim(only_predicted)[2]),
    "mse_Ross" = NA, "mse_East_Antarctica" = NA,
    "mse_Weddell" = NA, "mse_Bellingshausen_Amundsen_Sea" = NA,
    "mse_Total" = NA,
    # climatology
    "mse0_King_Haakon" = NA, "mse0_Ross" = NA, "mse0_East_Antarctica" = NA,
    "mse0_Weddell" = NA, "mse0_Bellingshausen_Amundsen_Sea" = NA,
    "mse0_Total" = NA)
  for (i in 1:dim(only_predicted)[2]) {
    res_data_err[i, 1:6] = (avg_prediction[i,] - this.all_dat[i,c(4:8, 3)])
    res_data_err[i,7:12] = (clim_avg_cv[i,] - this.all_dat[i,c(4:8, 3)])
  }
  cov_err = res_data_err[1:6] %>%
    cov()
  colnames(cov_err) = rownames(cov_err) = c("King H. VII", "Ross Sea", "East Ant.", "Weddell Sea", "B.A. Sea", "Total")
  round(cov_err * 10, 2)[c(2, 5, 4, 1, 3, 6),c(2, 5, 4, 1, 3, 6)]
  print(xtable::xtable(cov_err[c(2, 5, 4, 1, 3, 6),c(2, 5, 4, 1, 3, 6)] * 10))


  # covariance of sea ice anomaly
  cov_sie = cov(all_dat[,c(4:8,3)], use = "pair")
  colnames(cov_sie) = rownames(cov_sie) = c("King H. VII", "Ross Sea", "East Ant.", "Weddell Sea", "B.A. Sea", "Total")
  round(cov_sie * 10, 2)[c(2, 5, 4, 1, 3, 6),c(2, 5, 4, 1, 3, 6)]
  print(xtable::xtable(cov_sie[c(2, 5, 4, 1, 3, 6),c(2, 5, 4, 1, 3, 6)] * 10))

  # cor_sie = cor(all_dat[,c(4:8,3)], use = "pair")
  # round(cor_sie, 2)
  # plot sea ice anomaly
  png("../plots/sea_ice/sea_ice_ts.png", height = 360, width = 540)
  matplot(all_dat$tdate[958:1464], all_dat[958:1464,c(3,4:8,3)], type = "l", lty = 1,
          col = c("steelblue", rep("darkgrey", 5), "steelblue"), lwd = c(1.5,1,1,1,1,1,1.5),
          xlab = "Date",
          ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'))
  dev.off()

  my_rmse = function(x, y) {
    sqrt(mean((x - y)^2))
  }
  my_ae = function(x, y) {
    mean(abs(x - y))
  }

  # compile table of correlation, CE, and RMSE by sector for model and climatology
  res_table = data.frame(Sector = c("King H. VII", "Ross Sea", "East Ant.", "Weddell Sea", "B.A. Sea", "Total"),
                         cor_clim = c(cor(clim_avg_cv[,1], this.all_dat[,4]),
                                      cor(clim_avg_cv[,2], this.all_dat[,5]),
                                      cor(clim_avg_cv[,3], this.all_dat[,6]),
                                      cor(clim_avg_cv[,4], this.all_dat[,7]),
                                      cor(clim_avg_cv[,5], this.all_dat[,8]),
                                      cor(clim_avg_cv[,6], this.all_dat[,3])),
                         cor_mod = c(cor(avg_prediction[,1], this.all_dat[,4]),
                                     cor(avg_prediction[,2], this.all_dat[,5]),
                                     cor(avg_prediction[,3], this.all_dat[,6]),
                                     cor(avg_prediction[,4], this.all_dat[,7]),
                                     cor(avg_prediction[,5], this.all_dat[,8]),
                                     cor(avg_prediction[,6], this.all_dat[,3])),
                         ce_clim = c(my_CE(this.all_dat[,4], clim_avg_cv[,1]),
                                     my_CE(this.all_dat[,5], clim_avg_cv[,2]),
                                     my_CE(this.all_dat[,6], clim_avg_cv[,3]),
                                     my_CE(this.all_dat[,7], clim_avg_cv[,4]),
                                     my_CE(this.all_dat[,8], clim_avg_cv[,5]),
                                     my_CE(this.all_dat[,3], clim_avg_cv[,6])),
                         ce_mod = c(my_CE(this.all_dat[,4], avg_prediction[,1]),
                                    my_CE(this.all_dat[,5], avg_prediction[,2]),
                                    my_CE(this.all_dat[,6], avg_prediction[,3]),
                                    my_CE(this.all_dat[,7], avg_prediction[,4]),
                                    my_CE(this.all_dat[,8], avg_prediction[,5]),
                                    my_CE(this.all_dat[,3], avg_prediction[,6])),
                         ae_clim = c(my_ae(this.all_dat[,4], clim_avg_cv[,1]),
                                     my_ae(this.all_dat[,5], clim_avg_cv[,2]),
                                     my_ae(this.all_dat[,6], clim_avg_cv[,3]),
                                     my_ae(this.all_dat[,7], clim_avg_cv[,4]),
                                     my_ae(this.all_dat[,8], clim_avg_cv[,5]),
                                     my_ae(this.all_dat[,3], clim_avg_cv[,6])),
                         ae_mod = c(my_ae(this.all_dat[,4], avg_prediction[,1]),
                                    my_ae(this.all_dat[,5], avg_prediction[,2]),
                                    my_ae(this.all_dat[,6], avg_prediction[,3]),
                                    my_ae(this.all_dat[,7], avg_prediction[,4]),
                                    my_ae(this.all_dat[,8], avg_prediction[,5]),
                                    my_ae(this.all_dat[,3], avg_prediction[,6])),
                         rmse_clim = c(my_rmse(this.all_dat[,4], clim_avg_cv[,1]),
                                       my_rmse(this.all_dat[,5], clim_avg_cv[,2]),
                                       my_rmse(this.all_dat[,6], clim_avg_cv[,3]),
                                       my_rmse(this.all_dat[,7], clim_avg_cv[,4]),
                                       my_rmse(this.all_dat[,8], clim_avg_cv[,5]),
                                       my_rmse(this.all_dat[,3], clim_avg_cv[,6])),
                         rmse_mod = c(my_rmse(this.all_dat[,4], avg_prediction[,1]),
                                      my_rmse(this.all_dat[,5], avg_prediction[,2]),
                                      my_rmse(this.all_dat[,6], avg_prediction[,3]),
                                      my_rmse(this.all_dat[,7], avg_prediction[,4]),
                                      my_rmse(this.all_dat[,8], avg_prediction[,5]),
                                      my_rmse(this.all_dat[,3], avg_prediction[,6]))
  )#,
  # ae_clim = unlist(ae_dat[7:12]),
  # ae_mod = unlist(ae_dat[1:6]))
  res_table[,2:ncol(res_table)] = round(res_table[,2:ncol(res_table)],2)
  row.names(res_table) = NULL
  res_table[c(2, 5, 4, 1, 3, 6),]
  print(xtable::xtable(res_table[c(2, 5, 4, 1, 3, 6),]), include.rownames=FALSE)

  # compile table of correlation, CE, and RMSE by sector for model and climatology
  tmax = 420 # 2013
  # tmax = 384 # 2010
  
  res_table_2013 = data.frame(Sector = c("King H. VII", "Ross Sea", "East Ant.", "Weddell Sea", "B.A. Sea", "Total"),
                         cor_clim = c(cor(clim_avg_cv[1:tmax,1], this.all_dat[1:tmax,4]),
                                      cor(clim_avg_cv[1:tmax,2], this.all_dat[1:tmax,5]),
                                      cor(clim_avg_cv[1:tmax,3], this.all_dat[1:tmax,6]),
                                      cor(clim_avg_cv[1:tmax,4], this.all_dat[1:tmax,7]),
                                      cor(clim_avg_cv[1:tmax,5], this.all_dat[1:tmax,8]),
                                      cor(clim_avg_cv[1:tmax,6], this.all_dat[1:tmax,3])),
                         cor_mod = c(cor(avg_prediction[1:tmax,1], this.all_dat[1:tmax,4]),
                                     cor(avg_prediction[1:tmax,2], this.all_dat[1:tmax,5]),
                                     cor(avg_prediction[1:tmax,3], this.all_dat[1:tmax,6]),
                                     cor(avg_prediction[1:tmax,4], this.all_dat[1:tmax,7]),
                                     cor(avg_prediction[1:tmax,5], this.all_dat[1:tmax,8]),
                                     cor(avg_prediction[1:tmax,6], this.all_dat[1:tmax,3])),
                         ce_clim = c(my_CE(this.all_dat[1:tmax,4], clim_avg_cv[1:tmax,1]),
                                     my_CE(this.all_dat[1:tmax,5], clim_avg_cv[1:tmax,2]),
                                     my_CE(this.all_dat[1:tmax,6], clim_avg_cv[1:tmax,3]),
                                     my_CE(this.all_dat[1:tmax,7], clim_avg_cv[1:tmax,4]),
                                     my_CE(this.all_dat[1:tmax,8], clim_avg_cv[1:tmax,5]),
                                     my_CE(this.all_dat[1:tmax,3], clim_avg_cv[1:tmax,6])),
                         ce_mod = c(my_CE(this.all_dat[1:tmax,4], avg_prediction[1:tmax,1]),
                                    my_CE(this.all_dat[1:tmax,5], avg_prediction[1:tmax,2]),
                                    my_CE(this.all_dat[1:tmax,6], avg_prediction[1:tmax,3]),
                                    my_CE(this.all_dat[1:tmax,7], avg_prediction[1:tmax,4]),
                                    my_CE(this.all_dat[1:tmax,8], avg_prediction[1:tmax,5]),
                                    my_CE(this.all_dat[1:tmax,3], avg_prediction[1:tmax,6])),
                         ae_clim = c(my_ae(this.all_dat[1:tmax,4], clim_avg_cv[1:tmax,1]),
                                     my_ae(this.all_dat[1:tmax,5], clim_avg_cv[1:tmax,2]),
                                     my_ae(this.all_dat[1:tmax,6], clim_avg_cv[1:tmax,3]),
                                     my_ae(this.all_dat[1:tmax,7], clim_avg_cv[1:tmax,4]),
                                     my_ae(this.all_dat[1:tmax,8], clim_avg_cv[1:tmax,5]),
                                     my_ae(this.all_dat[1:tmax,3], clim_avg_cv[1:tmax,6])),
                         ae_mod = c(my_ae(this.all_dat[1:tmax,4], avg_prediction[1:tmax,1]),
                                    my_ae(this.all_dat[1:tmax,5], avg_prediction[1:tmax,2]),
                                    my_ae(this.all_dat[1:tmax,6], avg_prediction[1:tmax,3]),
                                    my_ae(this.all_dat[1:tmax,7], avg_prediction[1:tmax,4]),
                                    my_ae(this.all_dat[1:tmax,8], avg_prediction[1:tmax,5]),
                                    my_ae(this.all_dat[1:tmax,3], avg_prediction[1:tmax,6])),
                         rmse_clim = c(my_rmse(this.all_dat[1:tmax,4], clim_avg_cv[1:tmax,1]),
                                       my_rmse(this.all_dat[1:tmax,5], clim_avg_cv[1:tmax,2]),
                                       my_rmse(this.all_dat[1:tmax,6], clim_avg_cv[1:tmax,3]),
                                       my_rmse(this.all_dat[1:tmax,7], clim_avg_cv[1:tmax,4]),
                                       my_rmse(this.all_dat[1:tmax,8], clim_avg_cv[1:tmax,5]),
                                       my_rmse(this.all_dat[1:tmax,3], clim_avg_cv[1:tmax,6])),
                         rmse_mod = c(my_rmse(this.all_dat[1:tmax,4], avg_prediction[1:tmax,1]),
                                      my_rmse(this.all_dat[1:tmax,5], avg_prediction[1:tmax,2]),
                                      my_rmse(this.all_dat[1:tmax,6], avg_prediction[1:tmax,3]),
                                      my_rmse(this.all_dat[1:tmax,7], avg_prediction[1:tmax,4]),
                                      my_rmse(this.all_dat[1:tmax,8], avg_prediction[1:tmax,5]),
                                      my_rmse(this.all_dat[1:tmax,3], avg_prediction[1:tmax,6]))
  )#,
  # ae_clim = unlist(ae_dat[7:12]),
  # ae_mod = unlist(ae_dat[1:6]))
  res_table_2013[,2:ncol(res_table_2013)] = round(res_table_2013[,2:ncol(res_table_2013)],2)
  row.names(res_table_2013) = NULL
  res_table_2013[c(2, 5, 4, 1, 3, 6),]
  print(xtable::xtable(res_table_2013[c(2, 5, 4, 1, 3, 6),]), include.rownames=FALSE)


  sqrt(colMeans(all_dat[961:1464,c(4:8,3)]^2)) %>% round(2)

  # climatology (predict 0) vs model
  # climatology - model
  this.span = 0.025
  pdf(file = "../plots/BayesianModel/AE_climatology_model.pdf", width = 7, height = 5)
  plot(this.all_dat$tdate, res_data_abs[,7] - res_data_abs[,1], type = "l",
       lty = 1, main = "King Haakon VII", xlab = "Date", ylab = "AE climatology - AE model")
  abline(h = 0, lty = 2)
  lines(loess.smooth(this.all_dat$tdate, res_data_abs[,7] - res_data_abs[,1], span = this.span), col = 2)
  plot(this.all_dat$tdate, res_data_abs[,8] - res_data_abs[,2], type = "l",
       lty = 1, main = "Ross Sea", xlab = "Date", ylab = "AE climatology - AE model")
  abline(h = 0, lty = 2)
  lines(loess.smooth(this.all_dat$tdate, res_data_abs[,8] - res_data_abs[,2], span = this.span), col = 2)
  plot(this.all_dat$tdate, res_data_abs[,9] - res_data_abs[,3], type = "l",
       lty = 1, main = "East Antarctica", xlab = "Date", ylab = "AE climatology - AE model")
  abline(h = 0, lty = 2)
  lines(loess.smooth(this.all_dat$tdate, res_data_abs[,9] - res_data_abs[,3], span = this.span), col = 2)
  plot(this.all_dat$tdate, res_data_abs[,10] - res_data_abs[,4], type = "l",
       lty = 1, main = "Weddell Sea", xlab = "Date", ylab = "AE climatology - AE model")
  abline(h = 0, lty = 2)
  lines(loess.smooth(this.all_dat$tdate, res_data_abs[,10] - res_data_abs[,4], span = this.span), col = 2)
  plot(this.all_dat$tdate, res_data_abs[,11] - res_data_abs[,5], type = "l",
       lty = 1, main = "Bellingshausen Amundsen Sea", xlab = "Date", ylab = "AE climatology - AE model")
  abline(h = 0, lty = 2)
  lines(loess.smooth(this.all_dat$tdate, res_data_abs[,11] - res_data_abs[,5], span = this.span), col = 2)
  plot(this.all_dat$tdate, res_data_abs[,12] - res_data_abs[,6], type = "l",
       lty = 1, main = "Total", xlab = "Date", ylab = "AE climatology - AE model")
  abline(h = 0, lty = 2)
  lines(loess.smooth(this.all_dat$tdate, res_data_abs[,12] - res_data_abs[,6], span = this.span), col = 2)
  dev.off()

  #pdf(file = "../plots/performance/BayesianModel/AE_boxplots.pdf", width = 7, height = 5)
  boxplot(res_data_abs[,c(1,7,2,8,3,9,4,10,5,11,6,12)], col = 2:3,
          names = paste(rep(1:6, each = 2),c("Mod.", "Clim.")),
          at = c(1,2,4,5,7,8,10,11,13,14,16,17),
          ylim = c(0, 1), main = "Absolute error")
  legend("topleft", fill = 2:3, legend = c("Model", "Climatology"))
  dev.off()

  # rmse difference aggregated monthly, by year, 5 years, decade
  yearly_rmse = res_data2 %>%
    # monthly
    # mutate(year = 1:504) %>%
    # yearly
    # mutate(year = rep(1979:2020, each = 12)) %>%
    # 5 year blocks
    mutate(year = rep(1:6, each = 84)) %>%
    # decade
    # mutate(year = rep(c(rep(1, 11), rep(2:4, each = 10), 4), each = 12)) %>%
    group_by(year) %>%
    summarise_all(mean) %>%
    mutate_all(sqrt) %>%
    ungroup()
  # rmse difference
  my_diff = function(mod, clim) {
    (clim - mod)
  }
  yearly_rmse_diff = data.frame(
    Ross_Sea = my_diff(yearly_rmse$mse_Ross, yearly_rmse$mse0_Ross),
    Bellingsh_Amundsen = my_diff(yearly_rmse$mse_Bellingshausen_Amundsen_Sea,
                                 yearly_rmse$mse0_Bellingshausen_Amundsen_Sea),
    Weddell_Sea = my_diff(yearly_rmse$mse_Weddell, yearly_rmse$mse0_Weddell),
    King_Haakon = my_diff(yearly_rmse$mse_King_Haakon, yearly_rmse$mse0_King_Haakon),
    East_Antarctica = my_diff(yearly_rmse$mse_East_Antarctica, yearly_rmse$mse0_East_Antarctica),
    Total = my_diff(yearly_rmse$mse_Total, yearly_rmse$mse0_Total)
  )
  boxplot(yearly_rmse_diff, names = c("R.S.", "B.A.S.", "W.S", "K.H.S.", "E.A.", "Tot."))
  abline(h = 0, lty = 2)

  my_frac = function(mod, clim) {
    (clim - mod) / clim
  }
  # rmse fraction
  yearly_rmse_frac = data.frame(
    Ross_Sea = my_frac(yearly_rmse$mse_Ross, yearly_rmse$mse0_Ross),
    Bellingsh_Amundsen = my_frac(yearly_rmse$mse_Bellingshausen_Amundsen_Sea,
                                 yearly_rmse$mse0_Bellingshausen_Amundsen_Sea),
    Weddell_Sea = my_frac(yearly_rmse$mse_Weddell, yearly_rmse$mse0_Weddell),
    King_Haakon = my_frac(yearly_rmse$mse_King_Haakon, yearly_rmse$mse0_King_Haakon),
    East_Antarctica = my_frac(yearly_rmse$mse_East_Antarctica, yearly_rmse$mse0_East_Antarctica),
    Total = my_frac(yearly_rmse$mse_Total, yearly_rmse$mse0_Total)
  )
  boxplot(yearly_rmse_frac, names = c("R.S.", "B.A.S.", "W.S", "K.H.S.", "E.A.", "Tot."))
  abline(h = 0, lty = 2)



  # get monthly average mse
  res_data3 = res_data2
  res_data3$month = 1:12
  res_data3 = res_data3 %>%
    group_by(month) %>%
    summarize_all(list(mean = mean))
  # rmse by month
  pdf(file = "../plots/BayesianModel/RMSE_month.pdf", width = 7, height = 5)
  matplot(1:12, sqrt(res_data3[,c(2,8)]), type = "l", main = "King Haakon VII", xlab = "Month", ylab = expression("RMSE in 10^6 km^2"))
  matplot(1:12, sqrt(res_data3[,c(3, 9)]), type = "l", main = "Ross Sea", xlab = "Month", ylab = expression("RMSE in 10^6 km^2"))
  legend("topleft", fill = 1:2, legend = c("Model", "Climatology"))
  matplot(1:12, sqrt(res_data3[,c(4, 10)]), type = "l", main = "East Antarctica", xlab = "Month", ylab = expression("RMSE in 10^6 km^2"))
  matplot(1:12, sqrt(res_data3[,c(5, 11)]), type = "l", main = "Weddell Sea", xlab = "Month", ylab = expression("RMSE in 10^6 km^2"))
  matplot(1:12, sqrt(res_data3[,c(6, 12)]), type = "l", main = "Bellingshausen Amundsen Sea", xlab = "Month", ylab = expression("RMSE in 10^6 km^2"))
  matplot(1:12, sqrt(res_data3[,c(7, 13)]), type = "l", main = "Total", xlab = "Month", ylab = expression("RMSE in 10^6 km^2"))
  dev.off()

  # overall root mean squared error for model and climatology
  sqrt(mean(unlist(res_data2[,1:5])))
  sqrt(mean(unlist(res_data2[,6:10])))

  # rmse for total
  sqrt(mean(unlist(res_data2[,6])))
  sqrt(mean(unlist(res_data2[,12])))

  # # compute performance table for seasonal averages
  # # compute average prediction by sector and month
  # avg_prediction_seasonal = data.frame("King_Haakon" = NA, "Ross" = NA,
  #                                      "East_Antarctica" = NA, "Weddell" = NA,
  #                                      "Bellingshausen_Amundsen" = NA, "Total" = NA)
  # this.all_dat_seasonal = data.frame("King_Haakon" = NA, "Ross" = NA,
  #                                    "East_Antarctica" = NA, "Weddell" = NA,
  #                                    "Bellingshausen_Amundsen" = NA, "Total" = NA)
  # season = c(1,1, rep(2:168, each = 3), 169)
  # for (i in 1:169) {
  #   avg_prediction_seasonal[i, 1:6] = colMeans(avg_prediction[season == i,])
  #   this.all_dat_seasonal[i,1:6] = colMeans(this.all_dat[season == i,c(4:8,3)])
  # }
  # my_rmse = function(x, y) {
  #   sqrt(mean((x - y)^2))
  # }
  # my_ae = function(x, y) {
  #   mean(abs(x - y))
  # }
  #
  # # compile table of correlation, CE, and RMSE by sector for model and climatology
  # res_table_seasonal = data.frame(
  #   Sector = c("King Haakon VII", "Ross Sea", "East Antarctica", "Weddell Sea", "Amundsen Bellings.", "Total"),
  #   cor_clim = 0,
  #   cor_mod = c(cor(avg_prediction_seasonal[,1], this.all_dat_seasonal[,1]),
  #               cor(avg_prediction_seasonal[,2], this.all_dat_seasonal[,2]),
  #               cor(avg_prediction_seasonal[,3], this.all_dat_seasonal[,3]),
  #               cor(avg_prediction_seasonal[,4], this.all_dat_seasonal[,4]),
  #               cor(avg_prediction_seasonal[,5], this.all_dat_seasonal[,5]),
  #               cor(avg_prediction_seasonal[,6], this.all_dat_seasonal[,6])),
  #   ce_clim = 0,
  #   ce_mod = c(my_CE(this.all_dat_seasonal[,1], avg_prediction_seasonal[,1]),
  #              my_CE(this.all_dat_seasonal[,2], avg_prediction_seasonal[,2]),
  #              my_CE(this.all_dat_seasonal[,3], avg_prediction_seasonal[,3]),
  #              my_CE(this.all_dat_seasonal[,4], avg_prediction_seasonal[,4]),
  #              my_CE(this.all_dat_seasonal[,5], avg_prediction_seasonal[,5]),
  #              my_CE(this.all_dat_seasonal[,6], avg_prediction_seasonal[,6])),
  #   ae_clim =c(my_ae(this.all_dat_seasonal[,1], 0),
  #              my_ae(this.all_dat_seasonal[,2], 0),
  #              my_ae(this.all_dat_seasonal[,3], 0),
  #              my_ae(this.all_dat_seasonal[,4], 0),
  #              my_ae(this.all_dat_seasonal[,5], 0),
  #              my_ae(this.all_dat_seasonal[,6], 0)),
  #   ae_mod = c(my_ae(this.all_dat_seasonal[,1], avg_prediction_seasonal[,1]),
  #              my_ae(this.all_dat_seasonal[,2], avg_prediction_seasonal[,2]),
  #              my_ae(this.all_dat_seasonal[,3], avg_prediction_seasonal[,3]),
  #              my_ae(this.all_dat_seasonal[,4], avg_prediction_seasonal[,4]),
  #              my_ae(this.all_dat_seasonal[,5], avg_prediction_seasonal[,5]),
  #              my_ae(this.all_dat_seasonal[,6], avg_prediction_seasonal[,6])),
  #   rmse_clim =c(my_rmse(this.all_dat_seasonal[,1], 0),
  #                my_rmse(this.all_dat_seasonal[,2], 0),
  #                my_rmse(this.all_dat_seasonal[,3], 0),
  #                my_rmse(this.all_dat_seasonal[,4], 0),
  #                my_rmse(this.all_dat_seasonal[,5], 0),
  #                my_rmse(this.all_dat_seasonal[,6], 0)),
  #   rmse_mod = c(my_rmse(this.all_dat_seasonal[,1], avg_prediction_seasonal[,1]),
  #                my_rmse(this.all_dat_seasonal[,2], avg_prediction_seasonal[,2]),
  #                my_rmse(this.all_dat_seasonal[,3], avg_prediction_seasonal[,3]),
  #                my_rmse(this.all_dat_seasonal[,4], avg_prediction_seasonal[,4]),
  #                my_rmse(this.all_dat_seasonal[,5], avg_prediction_seasonal[,5]),
  #                my_rmse(this.all_dat_seasonal[,6], avg_prediction_seasonal[,6]))
  # )
  # res_table_seasonal[,2:ncol(res_table_seasonal)] = round(res_table_seasonal[,2:ncol(res_table_seasonal)],2)
  # row.names(res_table_seasonal) = NULL
  # res_table_seasonal[c(2, 5, 4, 1, 3, 6),]
  # print(xtable::xtable(res_table_seasonal[c(2, 5, 4, 1, 3, 6),]), include.rownames=FALSE)

}

# pdf(paste0("../output/f.",fname,".od.",vname,".pdf"))
# check_residuals(fit)
# #autoplot(forecast(object = fit,h = 12))
# dev.off()

# str(fit, 1)
library("dplyr")
fit$stanfit@sim %>% str(1)
# extract(fit$stanfit, "Y_latent_change") %>% str()

# plot one draw
# # latent data contains all missing data
# extract(fit$stanfit, "latent_data")[[1]][15, ,] %>%
#   matplot(type = "l")
# fit
rstan::extract(fit$stanfit, "fit")[[1]][15, ,] %>%
  matplot(type = "l")
rstan::extract(fit$stanfit, "Y_latent")[[1]][5, ,] %>%
  matplot(x = all_dat$tdate, type = "l")
abline(v = 1979, lty = 2)
legend("topright",
       legend = c("King Haakon", "Ross", "East Antarctica", "Weddell", "Bellingshausen Amundsen"),
       fill = 1:5)
# compute total sum
rstan::extract(fit$stanfit, "Y_latent")[[1]][5, ,] %>%
  rowSums() %>%
  plot(x = all_dat$tdate, type = "l")

# all_fit = rowSums(extract(fit$stanfit, "Y_latent")[[1]],
#                           dims = c(1,2))
# str(all_fit)

# get all fits for total sea ice
fits =  rstan::extract(fit$stanfit, "fit")[[1]]
# fits =  rstan::extract(fit$stanfit, "fit")[[1]]
# fits =  extract(fit$stanfit, "Y_latent")[[1]]
# ###############################################################
# # destandardize fits
# multiply by monthly standard deviation
for (i in 1:(dim(fits)[1])) {
  fits[i,,] = fits[i,,] * as.matrix(monthly_sd_df[,2:6])
}
# ###############################################################
# compute total fit
total_fits = matrix(NA, nrow = dim(fits)[1], ncol = dim(fits)[2])
for (i in 1:dim(fits)[1]) {
  total_fits[i,] = rowSums(fits[i,,])
}
str(total_fits)

# create data for github
reconstructions = array(NA, dim = c(2500, 1464, 6),
                        dimnames = list(paste0("reconstruction_", 1:2500), 
                                        paste0("t", all_dat$Year, "_", all_dat$Month),
                                        c("Total", "King_Hakon_VII",
                                          "Ross_Sea", "East_Antarctica",
                                          "Weddell_Sea", "Bellingshausen_Amundsen_Sea")))
reconstructions[,,2:6] = fits
reconstructions[,,1] = total_fits
dimnames(reconstructions)
# rowSums(reconstructions[2,,2:6]) - total_fits[2,]
# write.csv2(reconstructions, "../output/reconstructions.csv")
saveRDS(reconstructions, "../output/reconstructions.RDS")

write.csv(t(reconstructions[,,1]), "../output/reconstructions_total.csv")
write.csv(t(reconstructions[,,2]), "../output/reconstructions_king_haakon_VII.csv")
write.csv(t(reconstructions[,,3]), "../output/reconstructions_ross_sea.csv")
write.csv(t(reconstructions[,,4]), "../output/reconstructions_east_antarctica.csv")
write.csv(t(reconstructions[,,5]), "../output/reconstructions__weddell_sea.csv")
write.csv(t(reconstructions[,,6]), "../output/reconstructions_bellingshausen_amundsen_sea.csv")

# saveRDS(reconstructions[,,1], "../output/reconstructions_total.RDS")
# saveRDS(reconstructions[,,2], "../output/reconstructions_king_haakon_VII.RDS")
# saveRDS(reconstructions[,,3], "../output/reconstructions_ross_sea.RDS")
# saveRDS(reconstructions[,,4], "../output/reconstructions_east_antarctica.RDS")
# saveRDS(reconstructions[,,5], "../output/reconstructions__weddell_sea.RDS")
# saveRDS(reconstructions[,,6], "../output/reconstructions_bellingshausen_amundsen_sea.RDS")



# matplot(x = all_dat$tdate, y = t(total_fits), type = "l", lty = 1,
#         col = alpha(1, 0.3), main = "Total")
# matplot(x = all_dat$tdate, y = t(total_fits[c(T,rep(F,5)),]), type = "l", lty = 1,
#         col = alpha(1, 0.01), main = "Total")
# lines(all_dat$tdate, colMeans(total_fits), col = 2)

# # look at fit for first 2 years
# matplot(x = t(total_fits[c(T,rep(F,5)),1:24]), type = "l", lty = 1,
#         col = alpha(1, 0.05), main = "Total")
# lines(colMeans(total_fits)[1:24], col = 2)
# # January 1900
# abline(v = 13, lty = 3)
#
# mean(total_fits[,13])
# sd(total_fits[,13])
# summary(total_fits[,13])
#
# plot(apply(total_fits[,1:48], 2, sd))
# abline(v = 13, lty = 3)

plot(x = all_dat$tdate, y = total_fits[500,], type = "l", lty = 1, col = 1, main = "Total")
plot(x = all_dat$tdate, y = total_fits[1000,], type = "l", lty = 1, col = 1, main = "Total")
plot(x = all_dat$tdate, y = total_fits[1500,], type = "l", lty = 1, col = 1, main = "Total")
plot(x = all_dat$tdate, y = total_fits[2000,], type = "l", lty = 1, col = 1, main = "Total")

# sea ice by sector
# get all fits for total sea ice

# png(file = "../plots/BayesianModel/predictions%02d.png", width = 7*80, height = 5*80)
# King_Hakon
matplot(all_dat$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,1]), type = "l", lty = 1, col = alpha(1, 0.01),
        main = "King Haakon VII", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5),
        xlim = c(1899, 2021))
        # xlim = c(1977, 1980))
lines(all_dat$tdate[1:957], colMeans(fits[,1:957,1]), col = 2)
lines(all_dat$tdate, all_dat$King_Hakon, col = "steelblue")

# Ross
matplot(all_dat$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,2]), type = "l", lty = 1, col = alpha(1, 0.01),
        main = "Ross Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5),
        xlim = c(1899, 2021))
lines(all_dat$tdate[1:957], colMeans(fits[,1:957,2]), col = 2)
lines(all_dat$tdate, all_dat$Ross, col = "steelblue")

# East_Antarctica
matplot(all_dat$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,3]), type = "l", lty = 1, col = alpha(1, 0.01),
        main = "East Antarctica", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.75, 0.75),
        xlim = c(1899, 2021))
lines(all_dat$tdate[1:957], colMeans(fits[,1:957,3]), col = 2)
lines(all_dat$tdate, all_dat$East_Antarctica, col = "steelblue")

# Weddell
matplot(all_dat$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,4]), type = "l", lty = 1, col = alpha(1, 0.01),
        main = "Weddell Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5),
        xlim = c(1899, 2021))
lines(all_dat$tdate[1:957], colMeans(fits[,1:957,4]), col = 2)
lines(all_dat$tdate, all_dat$Weddell, col = "steelblue")

# Bellingshausen_Amundsen_Sea
matplot(all_dat$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,5]), type = "l", lty = 1, col = alpha(1, 0.01),
        main = "Bellinghausen Amundsen Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.75, 0.75),
        xlim = c(1899, 2021))
lines(all_dat$tdate[1:957], colMeans(fits[,1:957,5]), col = 2)
lines(all_dat$tdate, all_dat$Bellingshausen_Amundsen_Sea, col = "steelblue")

# Total
matplot(all_dat$tdate[1:957], t(total_fits[c(T,rep(F,5)),1:957]), type = "l", lty = 1, col = alpha(1, 0.01),
        main = "Total", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-2, 2),
        xlim = c(1899, 2021))
lines(all_dat$tdate[1:957], colMeans(total_fits)[1:957], col = 2)
# lines(all_dat$tdate[1:957], t(total_fits[1000,1:957]), type = "l", lty = 1, col = "royalblue")
lines(all_dat$tdate, all_dat$total, col = "steelblue")
dev.off()


# png(file = "../plots/BayesianModel/WeddellPresentation01.png", width = 7*80, height = 5*80)
# Weddell
plot(all_dat$tdate, all_dat$Weddell, type = "l", col = "steelblue",
     main = "Weddell Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
     ylim = c(-1.5, 1.5), xlim = c(1899, 2021))
dev.off()

# png(file = "../plots/BayesianModel/WeddellPresentation02.png", width = 7*80, height = 5*80)
# Weddell
plot(all_dat$tdate[1:957], t(fits[1000,1:957,4]), type = "l", lty = 1, col = "darkgrey", # alpha(1, 0.01),
        main = "Weddell Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5),
        xlim = c(1899, 2021))
# lines(all_dat$tdate[1:957], colMeans(fits[,1:957,4]), col = 2)
lines(all_dat$tdate, all_dat$Weddell, col = "steelblue")
dev.off()

# png(file = "../plots/BayesianModel/WeddellPresentation03.png", width = 7*80, height = 5*80)
# Weddell
matplot(all_dat$tdate[1:957], t(fits[c(T,rep(F,5)),1:957,4]), type = "l", lty = 1, col = alpha(1, 0.01),
     main = "Weddell Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
     ylim = c(-1.5, 1.5),
     xlim = c(1899, 2021))
# lines(all_dat$tdate[1:957], colMeans(fits[,1:957,4]), col = 2)
lines(all_dat$tdate, all_dat$Weddell, col = "steelblue")
dev.off()



# plot individual draws from posterior distribution
png(file = "../plots/BayesianModel/prediction_draws_%02d.png", width = 7*80, height = 5*80)
plot(all_dat$tdate[1:957], t(total_fits[c(100),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     main = "Total", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'), xlab = "Date")
lines(all_dat$tdate, all_dat$total, col = "steelblue")
abline(v = 1978.792, lty = 2)

plot(all_dat$tdate[1:957], t(total_fits[c(500),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     main = "Total", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'), xlab = "Date")
lines(all_dat$tdate, all_dat$total, col = "steelblue")
abline(v = 1978.792, lty = 2)

plot(all_dat$tdate[1:957], t(total_fits[c(900),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     main = "Total", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'), xlab = "Date")
lines(all_dat$tdate, all_dat$total, col = "steelblue")
abline(v = 1978.792, lty = 2)

plot(all_dat$tdate[1:957], t(total_fits[c(1200),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     main = "Total", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'), xlab = "Date")
lines(all_dat$tdate, all_dat$total, col = "steelblue")
abline(v = 1978.792, lty = 2)

plot(all_dat$tdate[1:957], t(total_fits[c(1600),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     main = "Total", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'), xlab = "Date")
lines(all_dat$tdate, all_dat$total, col = "steelblue")
abline(v = 1978.792, lty = 2)

plot(all_dat$tdate[1:957], t(total_fits[c(2000),1:957]), type = "l", lty = 1, col = c("steelblue"),
     xlim = c(1899, 2021), ylim = c(-2.2, 2.2),
     main = "Total", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'), xlab = "Date")
lines(all_dat$tdate, all_dat$total, col = "steelblue")
abline(v = 1978.792, lty = 2)
dev.off()


# plot all draws from posterior distribution for one year
png(file = "../plots/BayesianModel/prediction_1977.png", width = 7*80, height = 5*80)
matplot(all_dat$tdate[900:1000], t(total_fits[,900:1000]), type = "l", lty = 1, col = alpha(1, 0.05),
        xlim = c(1977.1, 1978.2), ylim = c(-1.75, 1.75),
        main = "Total", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'), xlab = "", xaxt="n")
axis(side = 1, tick = T,
     at = c(1977.042, 1977.125, 1977.208, 1977.292, 1977.375, 1977.458, 1977.542, 1977.625, 1977.708, 1977.792, 1977.875, 1977.958,
            1978.042, 1978.125, 1978.208, 1978.292),
     labels = c("Jan\n77", "Feb\n77", "Mar\n77", "Apr\n77", "May\n77", "Jun\n77", "Jul\n77", "Aug\n77", "Sep\n77", "Oct\n77", "Nov\n77", "Dec\n77",
                "Jan\n78", "Feb\n78", "Mar\n78", "Apr\n78"),
     padj = 0.5)
abline(h = 0, lty = 2, col = "lightblue", lwd = 2)
lines(all_dat$tdate[1:957], colMeans(total_fits)[1:957], col = 2, lwd = 2)
dev.off()
# ##############################################################################
# ######## plot seasonal reconstruction for JJA but do not show in paper
# ##############################################################################
JJA_dat = all_dat[all_dat$Month %in% c(6, 7, 8),]
JJA_dat = JJA_dat %>%
  group_by(Year) %>%
  summarise(total = mean(total, na.rm = TRUE),
            King_Hakon = mean(King_Hakon, na.rm = TRUE),
            Ross = mean(Ross, na.rm = TRUE),
            East_Antarctica = mean(East_Antarctica, na.rm = TRUE),
            Weddell = mean(Weddell, na.rm = TRUE),
            Bellingshausen_Amundsen_Sea = mean(Bellingshausen_Amundsen_Sea, na.rm = TRUE),
            tdate = mean(tdate))
JJA_fits_monthly = fits[,all_dat$Month %in% c(6, 7, 8),]
JJA_fits = JJA_fits_monthly[,1:122,]
for (i in 1:122) {
  for (j in 1:(dim(JJA_fits)[1])) {
    JJA_fits[j,i,] = apply(JJA_fits_monthly[j,(3 * i - 2):(3 * i),], 2, mean)
  }
}

# compute JJA total fit
JJA_total_fits = matrix(NA, nrow = dim(JJA_fits)[1], ncol = dim(JJA_fits)[2])
for (i in 1:dim(JJA_fits)[1]) {
  JJA_total_fits[i,] = rowSums(JJA_fits[i,,])
}
str(JJA_total_fits)

# png(file = "../plots/BayesianModel/predictionsJJA%02d.png", width = 7*80, height = 5*80)
# King_Hakon
matplot(JJA_dat$tdate[1:78], t(JJA_fits[c(T,rep(F,5)),1:78,1]), type = "l", lty = 1, col = alpha(1, 0.1),
        main = "King Haakon VII", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5),
        xlim = c(1899, 1980))
# xlim = c(1977, 1980))
# lines(all_dat$tdate[1:957], colMeans(fits[,1:957,1]), col = 3)
lines(JJA_dat$tdate[1:78], colMeans(JJA_fits[,1:78,1]), col = 2, lwd = 2)

# Ross
matplot(JJA_dat$tdate[1:78], t(JJA_fits[c(T,rep(F,5)),1:78,2]), type = "l", lty = 1, col = alpha(1, 0.1),
        main = "Ross Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5),
        xlim = c(1899, 1980))
# lines(all_dat$tdate[1:957], colMeans(fits[,1:957,2]), col = 3)
lines(JJA_dat$tdate[1:78], colMeans(JJA_fits[,1:78,2]), col = 2, lwd = 2)

# East_Antarctica
matplot(JJA_dat$tdate[1:78], t(JJA_fits[c(T,rep(F,5)),1:78,3]), type = "l", lty = 1, col = alpha(1, 0.1),
        main = "East Antarctica", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.75, 0.75),
        xlim = c(1899, 1980))
# lines(all_dat$tdate[1:957], colMeans(fits[,1:957,3]), col = 3)
lines(JJA_dat$tdate[1:78], colMeans(JJA_fits[,1:78,3]), col = 2, lwd = 2)

# Weddell
matplot(JJA_dat$tdate[1:78], t(JJA_fits[c(T,rep(F,5)),1:78,4]), type = "l", lty = 1, col = alpha(1, 0.1),
        main = "Weddell Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-1.5, 1.5),
        xlim = c(1899, 1980))
# lines(all_dat$tdate[1:957], colMeans(fits[,1:957,4]), col = 3)
lines(JJA_dat$tdate[1:78], colMeans(JJA_fits[,1:78,4]), col = 2, lwd = 2)

# Bellingshausen_Amundsen_Sea
matplot(JJA_dat$tdate[1:78], t(JJA_fits[c(T,rep(F,5)),1:78,5]), type = "l", lty = 1, col = alpha(1, 0.1),
        main = "Bellinghausen Amundsen Sea", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-0.75, 0.75),
        xlim = c(1899, 1980))
# lines(all_dat$tdate[1:957], colMeans(fits[,1:957,5]), col = 3)
lines(JJA_dat$tdate[1:78], colMeans(JJA_fits[,1:78,5]), col = 2, lwd = 2)

# Total
matplot(JJA_dat$tdate[1:78], t(JJA_total_fits[c(T,rep(F,5)),1:78]), type = "l", lty = 1, col = alpha(1, 0.1),
        main = "Total", xlab = "Date", ylab = parse(text='"Sea ice anomaly in 10^6 km^2"'),
        ylim = c(-2, 2),
        xlim = c(1899, 1980))
# lines(all_dat$tdate[1:957], colMeans(total_fits)[1:957], col = 3)
lines(JJA_dat$tdate[1:78], colMeans(JJA_total_fits)[1:78], col = 2, lwd = 2)
# lines(all_dat$tdate[1:957], t(total_fits[1000,1:957]), type = "l", lty = 1, col = "royalblue")
dev.off()


# ################################################################################
# look at parameter summary
bayesbackcast:::summary.varstan(fit, pars = c("beta0C", "brg"))
# reformat to 5x5
bayesbackcast:::summary.varstan(fit, pars = "Sigma")
bayesbackcast:::summary.varstan(fit, pars = "Sigma") %>%
  select(mean) %>%
  unlist() %>%
  matrix(ncol = 5) %>%
  cov2cor() %>%
  round(2)

eps = bayesbackcast:::summary.varstan(fit, pars = "epsilon")
head(eps)
tail(eps)
eps %>%
  mutate(sector = rep(1:5, each = nrow(eps) / 5)) %>%
  group_by(sector) %>%
  summarise_all(mean, na.rm = TRUE)
str(eps)

extract(fit$stanfit, "epsilon")[[1]][100,,] %>%
  cor() %>% round(2)
extract(fit$stanfit, "Sigma")[[1]][100,,] %>%
  unlist() %>%
  matrix(ncol = 5) %>%
  cov2cor() %>%
  round(2)

# source("plot.roots.R")
# pdf(paste0("mv_",fname,".od.",vname,".pdf"))
# c01 <- summary(fit)
# c01_arima <- c01[c(grep('ar.1',rownames(c01),fixed=T),grep('ma.1',rownames(c01),fixed=T)),]
# c01_arima
# ar_roots <- c01_arima[1:Kar,1]
# ar_roots=polyroot(c(1,-ar_roots))
# #ma_roots <- c01_arima[Kar+(1:Kma),1]
# #ma_roots=polyroot(c(1,ma_roots))
# plot.armaroots(ar_roots,type="AR")
# #plot.armaroots(ma_roots,type="MA")
# abs(ar_roots) # stationary iff all greater than 1
# #abs(1/ar_roots) # stationary iff all less than 1
# #abs(ma_roots) # invertible iff all greater than 1

c01_beta <- bayesbackcast:::summary.varstan(fit, pars = c("beta0C", "brg")) #  summary(fit, pars=colnames(fit$model$X))
c01_beta[c01_beta[,3] > 0 | c01_beta[,4] < 0,]
c01_beta[order(-abs(c01_beta[,1]))[1:100],]

beta0C = extract(fit$stanfit, c("beta0C"))[[1]]
str(beta0C)
head(beta0C)
hist(beta0C[,1])
hist(beta0C[,2])
hist(beta0C[,3])
hist(beta0C[,4])
hist(beta0C[,5])

# n samples x n sea ice dims x n coefficients
brg = rstan::extract(fit$stanfit, c("brg"))[[1]]
hist(brg[,1,10])
sapply(brg[,1,12], all.equal, current = 0, tolerance = 0.01) %>% table()
sapply(as.vector(brg), all.equal, current = 0, tolerance = 0.01) %>% table()
# get coefficients for King Haakon sector
brg_KingHaakon = brg[,1,]
# for (i in 1:50) {
#   hist(brg_KingHaakon[,i], main = i, xlim = c(-0.1, 0.1))
# }
# # 70 and 87 are notable for diff = 1
# (apply(brg_KingHaakon, 2, sd) * 1000) %>% round()
# hist(brg_KingHaakon[,70])
# hist(brg_KingHaakon[,87])
# # it's either PC 70 or PC 87 that seems to be relevant
# plot(brg_KingHaakon[,70], brg_KingHaakon[,87])
#
# all_pc[,70] %>% plot(type = "l")
# abline(h = 0, col = 2)
# lm(all_pc[,70] ~ I(1:1464))
#
# all_pc[,87] %>% plot(type = "l")
# abline(h = 0, col = 2)
# abline(lm(all_pc[,87] ~ I(1:1464)), col = 2)
# lm(all_pc[,87] ~ I(1:1464))

# # get coefficients for Bellingshausen Amundsen sector
# brg_BellingshausenAmundsen = brg[,5,]
# colnames(brg_BellingshausenAmundsen) = names(all_x[,1:984])
# cormat = round(cor(abs(brg_BellingshausenAmundsen)), 2)
# # 890500 Bellingshausen
# colnames(cormat)[c(75:76, 535:544)]
# cormat[order(abs(cormat[,75])),75]
# cormat[order(abs(cormat[,76])),76]
# # 889030 Grytviken
# colnames(cormat)[c(65:66, 485:494)]
# cormat[c(75:76, 535:544, 65:66, 485:494), c(75:76, 535:544, 65:66, 485:494)] %>%
#   round(2) %>%
#   image()

# ##############################################################################
# compute some kind of variable importance
# avg_prediction = data.frame("King_Haakon" = NA, "Ross" = NA,
#                             "East_Antarctica" = NA, "Weddell" = NA,
#                             "Bellingshausen_Amundsen" = NA)
str(brg)
var_imp = matrix(NA, nrow = dim(brg)[2], ncol = dim(brg)[3])
for (i in 1:(dim(brg)[2])) {
  var_imp[i, ] = colMeans(abs(brg)[,i,])
}
colnames(var_imp) = names(all_x[,1:(ncol(all_x)-1)])
# colnames(var_imp) = substr(names(all_x[,1:984]), 4, 100)
# recombine all lags with original variable
var_imp_total = var_imp[,1:164]
match_ID = sapply(colnames(var_imp_total), function(x) {
  strsplit(x, split = "_l")[[1]][1]
})
for (i in 1:164) {
  var_imp_total[,i] = rowSums(var_imp[,match_ID == colnames(var_imp)[i]])
}
var_imp_total

# which(!(colnames(var_imp_total) %in% all_stations$ID))
# colnames(var_imp_total)[!(colnames(var_imp_total) %in% all_stations$ID)]

# read in station data
all_stations = read.csv("../data/fogt_predictor_data/station_data.csv")
# subset to stations that are in the model
all_stations = all_stations[all_stations$ID %in% colnames(var_imp_total),]

all_stations$var = sapply(all_stations$ID, function(x) {
  strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])]
})
# add fake stations for climate indices
index_stations = data.frame(ID = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                   "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                   "SSTsNino4Mean.1854"),
                            longname = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                         "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                         "SSTsNino4Mean.1854"),
                            name = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                     "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                     "SSTsNino4Mean.1854"),
                            lat = -85, lon = seq(-150,150, length.out = 9),
                            ID.1 = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                     "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                     "SSTsNino4Mean.1854"),
                            start = "1900", end = "2019", var = "ind")
all_stations = rbind(all_stations, index_stations)

# check that all names are alphabetically ordered
order(colnames(var_imp_total))
order(all_stations$ID.1)

# add importance
station_importance = all_stations
station_importance[,c("importance_King_Haakon", "importance_Ross",
                      "importance_East_Antarctica", "importance_Weddell",
                      "importance_Bellingshausen_Amundsen")] = t(var_imp_total)

# aggregate by station
station_importance_agg = station_importance %>%
  group_by(name, lat, lon) %>%
  summarize(importance_King_Haakon = sum(importance_King_Haakon),
            importance_Ross = sum(importance_Ross),
            importance_East_Antarctica = sum(importance_East_Antarctica),
            importance_Weddell = sum(importance_Weddell),
            importance_Bellingshausen_Amundsen = sum(importance_Bellingshausen_Amundsen))

# plot station importance
ant <- ggplot2::map_data(map = "world", region = ".",
                         # orientation = c(-90, 0, 0),
                         wrap = c(-180, 180, -90))
sectors = c("King_Haakon", "Ross",
            "East_Antarctica", "Weddell",
            "Bellingshausen_Amundsen")
pretty_sectors = c("King Haakon VII", "Ross Sea",
                   "East Antarctica", "Weddell Sea",
                   "Bellingshausen Amundsen Sea")
for (i in 1:5) {
  this.sector = sectors[i]
  this.pretty_sector = pretty_sectors[i]
  # make sure we get the right importance per sector
  filter_dat = station_importance_agg[,c(1:3, 3+i)]
  names(filter_dat)[4] = "importance"
  filter_dat = filter_dat %>% ungroup() %>% filter(importance > quantile(importance, 0.75))

  density_dat = density(station_importance_agg$lon,
                        bw = 10,
                        from = -180, to = 180,
                        weights = unlist(station_importance_agg[,3+i]) / sum(unlist(station_importance_agg[,3+i])))
  density_dat = density(filter_dat$lon,
                        bw = 10,
                        from = -180, to = 180,
                        weights = unlist(filter_dat[,4]) / sum(unlist(filter_dat[,4])))
  density_dat = data.frame(x = density_dat$x, y = density_dat$y)
  # geographic projection
  my_col = rep("steelblue", 5)
  my_col[i] = "firebrick"
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group),
                 data = ant,
                 fill = "gray70", colour = "gray70") +
    # coord_map("stereographic", orientation=c(-90, 0, 0), ylim = c(-90,-50)) +
    scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(name = "latitude") +
    coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 0)) +
    # sea ice based sectors
    geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
    # sea ice based sectors
    geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5)[c(5,2,1,4,3)], y = -3,
                                 label = c("East\nAntarctica",
                                           "Ross Sea",
                                           "Amund.\nBel. Sea",
                                           "Weddell\nSea",
                                           "King Haakon")[c(5,2,1,4,3)]),
               aes(x = x, y = y, label = label),
               color = my_col) +
    theme_minimal() +
    # add stations
    geom_point(aes(x = lon, y = lat, size = importance, alpha = importance),
               data = filter_dat, show.legend = FALSE) +
    ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = importance),
                             max.overlaps = 20,
                             show.legend = FALSE,
                             data = filter_dat) +
    ggtitle(paste0(this.pretty_sector, ", ", "25% most important stations")) +
    scale_alpha("Variable\nImportance") +
    scale_size("Variable\nImportance") +
    geom_line(aes(x = x, y = y * 2000 - 85), data = density_dat)
  ggsave(filename = paste0("../plots/BayesianModel/StationImp_",
                           this.sector, ".png"),
         width = 8, height = 5,
         device = png())
  dev.off()

  filter_dat2 = slice_sample(filter_dat, n = 200, weight_by = filter_dat$importance, replace = TRUE)
  table(filter_dat2$name)

  ggplot() +
    # background density
    stat_density_2d(aes(x=lon, y=lat, fill = ..density..),
                    data = filter_dat2, geom = "raster", contour = FALSE, show.legend = FALSE) +
    # scale_fill_distiller(palette="Spectral", direction=-1) +
    scale_fill_distiller(palette="RdBu", direction=-1) +
    # scale_fill_distiller(palette="Blues", direction=1) +
    geom_polygon(aes(x = long, y = lat, group = group),
                 data = ant,
                 fill = NA,
                 colour = "gray70") +
    scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(name = "latitude") +
    coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 0)) +
    # sea ice based sectors
    geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
    # sea ice based sectors
    geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5)[c(5,2,1,4,3)], y = -3,
                                 label = c("East\nAntarctica",
                                           "Ross Sea",
                                           "Amund.\nBel. Sea",
                                           "Weddell\nSea",
                                           "King Haakon")[c(5,2,1,4,3)]),
               aes(x = x, y = y, label = label),
               color = my_col) +
    theme_minimal() +
    # add stations
    geom_point(aes(x = lon, y = lat, size = importance, alpha = importance),
               data = filter_dat, show.legend = FALSE) +
    ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = importance),
                             max.overlaps = 20,
                             show.legend = FALSE,
                             data = filter_dat) +
    ggtitle(paste0(this.pretty_sector, ", ", "25% most important stations")) +
    scale_alpha("Variable\nImportance") +
    scale_size("Variable\nImportance") +
    geom_line(aes(x = x, y = y * 2000 - 85), data = density_dat)
  ggsave(filename = paste0("../plots/BayesianModel/StationImpFill_",
                           this.sector, ".png"),
         width = 8, height = 5,
         device = png())
  dev.off()
}


# nu
nu = extract(fit$stanfit, c("nu"))[[1]]
# check nu for King Haakon
matplot(nu[c(T,F,F,F,F),,1],
        type = "l", lty = 1, col = alpha(1, 0.01))
matplot(t(nu[1:4,,1]),
        type = "l", lty = 1)
matplot(t(nu[1000:1003,,1]),
        type = "l", lty = 1)
matplot(t(nu[1996:2000,,1]),
        type = "l", lty = 1)

ar = extract(fit$stanfit, c("ar"))[[1]]
dim(ar)
# for King Haakon sector
hist(ar[,1,1])
hist(ar[,1,2])
hist(ar[,1,3])
hist(ar[,1,4])

# average ar coefficients per sector
colMeans(ar[,1,])
boxplot(ar[,1,])

colMeans(ar[,2,])
boxplot(ar[,2,])

colMeans(ar[,3,])
boxplot(ar[,3,])

colMeans(ar[,4,])
boxplot(ar[,4,])

colMeans(ar[,5,])
boxplot(ar[,5,])

boxplot(ar[,,1])
boxplot(ar[,,2])
boxplot(ar[,,3])
boxplot(ar[,,4])


# extract seasonal ar coefficients
ar_seasonal = extract(fit$stanfit, c("ar_seasonal"))[[1]]
dim(ar_seasonal)

# 1st sector lag 1 AR coefficient over the year
colMeans(ar_seasonal[,,1,1]) %>% plot()
boxplot(ar_seasonal[,,1,1])

# check if it is sensitive to my prior specification
# Y_t = AR1 * y_{t-1} + \varepsilon
# -1 < AR1 < 1

# 2nd sector lag 1 AR coefficient over the year
colMeans(ar_seasonal[,,2,1]) %>% plot()
boxplot(ar_seasonal[,,2,1])
# 3rd sector lag 1 AR coefficient over the year
colMeans(ar_seasonal[,,3,1]) %>% plot()
boxplot(ar_seasonal[,,3,1])
# 4th sector lag 1 AR coefficient over the year
colMeans(ar_seasonal[,,4,1]) %>% plot()
boxplot(ar_seasonal[,,4,1])
# 5th sector lag 1 AR coefficient over the year
colMeans(ar_seasonal[,,5,1]) %>% plot()
boxplot(ar_seasonal[,,5,1])

dim(ar_seasonal)
(cov(ar_seasonal[,,3,1])*100) %>% round(2) %>% image()


MASS::mvrnorm(n = 1000, mu = rep(0, 3),
              Sigma = matrix(c(.1, .05, 0,
                               .05, .1, .05,
                               0, .05, .1),
                             nrow = 3)) %>% cor() %>% round(2)



# average ma coefficients per sector
ma = extract(fit$stanfit, c("ma"))[[1]]
colMeans(ma)

# pars1 <- c("ar[1,1]","ma[1,1]","ar[1,2]","ma[1,2]","sigma0","beta_0")
# pars1 <- c("ar[1,1]","ma[1,1]","sigma0","beta_0")
# pars1 <- c("ar[1,1]","sigma0","beta_0")
check_hmc_diagnostics(fit$stanfit)
color_scheme_set("red")
#mcmc_plot(fit, pars = "sigma0")
#mcmc_pairs(fit$stanfit, pars = c("sigma0","beta_0"))
#mcmc_pairs(fit$stanfit, pars = pars1)
#mcmc_pairs(fit$stanfit, pars=colnames(fit$model$X)[1:5])
# colnames(fit$model$X)[1:5]
# samp <- as.matrix(extract(fit$stanfit, pars = colnames(fit$model$X)))
# samp[[1]][1:20]
#samp <- as.matrix(extract(fit$stanfit, pars = "brg"))
#samp[[1]][1:20]
# a <- sapply(samp,function(x){mean(abs(x))})
# a <- sapply(samp,function(x){median(x[abs(x)>0.0001])})
# c01_beta[order(-a)[1:20],]
# b <- rownames(samp)[order(-a)]
# mcmc_pairs(fit$stanfit, pars=b[1:5])
#g <- as.numeric(samp[(order(-a))[1],])[[1]]
# plot(density(as.numeric(samp[(order(-a))[1],][[1]])))
# plot(density(as.numeric(samp[(order(-a))[2],][[1]])))
# plot(density(as.numeric(samp[(order(-a))[3],][[1]])))
# plot(density(as.numeric(samp[(order(-a))[4],][[1]])))
# plot(density(as.numeric(samp[(order(-a))[5],][[1]])))
#print(fit, pars = pars, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
##pairs(stan_fit_sie pars = pars)
##traceplot(stan_fit_sie inc_warmup = TRUE, pars = pars)
###plot_draws(stan_fit_sie n_sample = 20, data_count = data.frame(X=X,y=Y_asy))
##plot(stan_fit_sie pars = pars)

# check 216 notes
# degree of model misspecification
if(F){
  # loo1 <- loo(fit, cores = 8)
  # pdf(paste0("../output/g.",fname,".od.",vname,".pdf"))
  # plot(loo1)
  # dev.off()
  # print(loo1)
  # str(loo1)

  # this is the loo on the observed years only
  LLarray <- loo::extract_log_lik(stanfit = fit$stanfit,
                                  parameter_name = "log_lik",
                                  merge_chains = FALSE)
  LLarray = LLarray[,,961:1464]
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = 10)
  this.loo <- loo::loo.array(LLarray, r_eff = r_eff, cores = 10,
                             save_psis = FALSE)
  plot(this.loo)
  this.loo
}
# # p_nonzero = 15, slab_scale = 0.025 with seasonally varying AR coefficient
# Computed from 2500 by 504 log-likelihood matrix
#
# Estimate   SE
# elpd_loo  -2070.7 46.0
# p_loo       271.6 10.4
# looic      4141.4 92.1
# ------
#   Monte Carlo SE of elpd_loo is NA.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     469   93.1%   6
# (0.5, 0.7]   (ok)        30    6.0%   1
# (0.7, 1]   (bad)        4    0.8%   0
# (1, Inf)   (very bad)   1    0.2%   0
# See help('pareto-k-diagnostic') for details.

# # p_nonzero = 15, slab_scale = 0.025 with seasonally varying AR coefficient
# elpd is less negative, so better
# Computed from 2500 by 504 log-likelihood matrix
#
# Estimate   SE
# elpd_loo  -2023.7 44.5
# p_loo       427.5 14.5
# looic      4047.4 88.9
# ------
#   Monte Carlo SE of elpd_loo is NA.
#
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     348   69.0%   6
# (0.5, 0.7]   (ok)       123   24.4%   1
# (0.7, 1]   (bad)       31    6.2%   1
# (1, Inf)   (very bad)   2    0.4%   4
# See help('pareto-k-diagnostic') for details.


# if(any(is.na(all_dat$sie))){
#   y_fit = extract(fit$stanfit, "fit")[[1]][15, ,]
#   matplot(all_dat$tdate, y_fit, type = "b",
#           xlim = c(1974, 1980), ylim = c(-1.5, 1.5))
#   y_latent = extract(fit$stanfit, "Y_latent")[[1]][15, ,]
#   matplot(all_dat$tdate, y_latent, type = "b",
#           xlim = c(1974, 1980), ylim = c(-1.5, 1.5))
#
#   residuals = extract(fit$stanfit, "residuals")[[1]][15, ,]
#   y_fit - (y_latent - residuals)
#
# }

dev.off()

# ##############################################################################
# how unusual was the 2014-2017 decline in a historic context?
plot(all_dat$tdate[all_dat$Year>2010], all_dat$total[all_dat$Year>2010], type = "l")
# the largest anomaly was in January 2015 at 1.943735 the smallest in December 2016 at -2.1844603942
# check how often a jump of 4.128195 or more occured within 2 years
total_fits %>% str()
which(all_dat$Year==1979 & all_dat$Month == 12)
max_jump = apply(total_fits[,1:972], 1, function(x) {
  this.max = 0
  for (i in 1:length(x)) {
    this.window = min(length(x), i + 36) # 24
    if (abs(max(x[i:this.window]) - min(x[i:this.window])) > this.max) {
      this.max = abs(max(x[i:this.window]) - min(x[i:this.window]))
    }
  }
  return(this.max)
})
str(max_jump)
# png("../plots/BayesianModel/p_value_decline.png", height = 400, width = 600)
hist(max_jump, xlim = c(0.5, 5), xlab = "3 year maximal difference", main = "")
abline(v = 4.128195, lty = 2)
dev.off()
mean(max_jump > 4.128195)
sum(max_jump > 4.128195)

total_fits %>% str()

# analyze the February 2022 sea ice minmum #####################################
# plot all the February estimates
observed_february = total_fits[,973:1464][1,all_dat$Month[973:1464] == 2]
hist(observed_february)
summary(observed_february)
# timeDate::skewness(observed_february)

# all_february = total_fits[,1:972][,all_dat$Month[1:972] == 2]
all_february = total_fits[,13:960][,all_dat$Month[13:960] == 2]
summary(all_february)
dim(all_february)

febyr = all_dat[13:960, "Year"][all_dat$Month[13:960] == 2]

# what's up with this year?
# all_february = all_february[,-78]
# get the minimum February estimates
all_february_minDraw = apply(all_february, 1, min)
# all_february_minDraw = apply(all_february, 1, quantile, 0.025)
# all_february_minDraw = apply(all_february, 1, quantile, 0.01)
# plot(apply(all_february, 2, sd))

# a good chunk of the reconstructions think that if a record sea ice minimum was achieved that it must have happened in
# 1977? for some reason
apply(all_february, 1, which.min) %>% table()
whichmin = apply(all_february, 1, which.min) + 1899
table(whichmin)
mean(whichmin == 1977)

png("../plots/BayesianModel/lowest_february.png", height = 400, width = 600)
hist(whichmin, breaks = seq(1900, 1980, 5),
     main = "", xlab = "Year of lowest reconstructed February anomaly")
dev.off()
# the anomaly is -0.9086785543
table(whichmin[all_february_minDraw < -0.9086785543])
# 105 out of 400 (26.25%) of draws where the minimum was below -0.96 were in 1977
mean(whichmin[all_february_minDraw < -0.9086785543] == 1977)
sum(whichmin[all_february_minDraw < -0.9086785543] == 1977)
hist(whichmin[all_february_minDraw < -0.9086785543],
     breaks = seq(1900, 1980, 5),
     xlab = "Year", main = "")


# png("../plots/BayesianModel/p_value_february.png", height = 400, width = 600)
hist(all_february_minDraw, main = "", xlab = "Lowest February Anomaly", xlim = c(-1.5, 0))
abline(v = -0.9086785543, lty = 2, lwd = 2)
abline(v = -1.088333, lty = 3, lwd = 2)
dev.off()
# 16% 
png("../plots/BayesianModel/p_value_february2.png", height = 400, width = 600)
hist(all_february_minDraw, main = "", xlab = "Lowest February Anomaly", xlim = c(-1.5, 0))
abline(v = -0.9086785543, lty = 2, lwd = 2)
abline(v = -1.088333, lty = 3, lwd = 2)
text(x = c(-0.9086785543 + 0.055, -1.088333 + 0.055), y = c(500, 500), 
     labels = c("2022", "2023"))
dev.off()
mean(all_february_minDraw < -0.9086785543)
sum(all_february_minDraw < -0.9086785543)

# february 2023
mean(all_february_minDraw < -1.088333)
sum(all_february_minDraw < -1.088333)
sum(whichmin[all_february_minDraw < -1.088333] == 1977)
mean(whichmin[all_february_minDraw < -1.088333] == 1977)


# check if march needs to be included
all_march = total_fits[,13:960][,all_dat$Month[13:960] == 3]
apply(all_march, 1, which.min) %>% table()
whichmin_march = apply(all_march, 1, which.min) + 1899
table(whichmin_march)

all_march_minDraw = apply(all_march, 1, min) - 3.12 + 4.14
hist(all_march_minDraw)
mean(all_march_minDraw < -0.96)
sum(all_march_minDraw < -0.96)
which(all_march_minDraw < -0.96)

# ##############################################################################
# // prior estimtaed of total sea ice extent
# // target += normal_lpdf(september1964|0.9,0.9); // Meier et al 2013 19.7 mioskm, range 18.8, 20.4, sd=0.8
# // target += normal_lpdf(may1966|0.4,0.2); // Gallager et al 2014 10.7 mioskm, sd 0.2
# // target += normal_lpdf(june1966|1,0.4); // 14.3 mioskm, sd 0.4
# // target += normal_lpdf(july1966|0.6,0.5); // 16.4 mioskm, sd 0.4
# // target += normal_lpdf(august1966|-1.6,0.4); // 15.9 mioskm,sd 0.4 // this august value is doubted

# I do not believe the FNN::KL.divergence KL divergence function to be correct
# I write my own, assuming both of them are normal
my_KL = function(x, y) {
  mu_x = mean(x)
  mu_y = mean(y)
  sd_x = sd(x)
  sd_y = sd(y)
  log(sd_y / sd_x) + (sd_x^2 + (mu_x - mu_y)^2) / (2 * sd_y^2) - 0.5
}
my_KL_test = replicate(10000, my_KL(rnorm(10000), rnorm(200)))
# hist(2*10 * my_KL_test)
quantile(2*200 * my_KL_test, probs = c(0.9, 0.95))
qchisq(c(0.9, 0.95), df = 2)

# for fixed mu_y and sd_y
my_KL2 = function(x, mu_y, sd_y) {
  mu_x = mean(x)
  sd_x = sd(x)
  log(sd_y / sd_x) + (sd_x^2 + (mu_x - mu_y)^2) / (2 * sd_y^2) - 0.5
}
my_KL_test = replicate(10000, my_KL2(x = rnorm(200), mu_y = 0, sd_y = 1))
quantile(2*200 * my_KL_test, probs = c(0.9, 0.95))
qchisq(c(0.9, 0.95), df = 2)

my_KL_test = replicate(10000, my_KL2(x = rnorm(1500), mu_y = 0, sd_y = 1))
quantile(2*1500 * my_KL_test, probs = c(0.9, 0.95))
qchisq(c(0.9, 0.95), df = 2)

# # perform a likelihood ratio test
# lik1 = pnorm()
# lik2 = pnorm()
# lambda = lik1 / lik2
# # compare to
# exp(-0.5*qchisq(0.95, df = 2))


# what does our model say in comparison
# should we use the kullback leibler divergence between their estimate and our posterior?
# should we just give a p value for our average reconstruction?
# september1964
# KLseptember1964 = FNN::KL.divergence(total_fits[,789], rnorm(n = 10000, mean = 0.9, sd = 0.4))[5]
KLseptember1964 = my_KL2(x = total_fits[,789], mu_y = 0.9, sd_y = 0.4)
2*1*KLseptember1964 # I pretend my sample size of the reconstruction is only one
qchisq(c(0.9, 0.95), df = 1)
qchisq(c(0.9, 0.95), df = 2)
# this is the p value
pchisq(2*1*KLseptember1964, df = 2)

# hist(rnorm(n = 500, mean = 0.9, sd = 0.9))
# abline(v=mean(total_fits[,789]))
dens = data.frame(x = seq(-3, 3, length.out = 100),
                  y = dnorm(seq(-3, 3, length.out = 100), mean = 0.9, sd = 0.4))
ggplot() +
  geom_density(aes(x = total_fits[,789])) +
  geom_line(aes(x = x, y = y), data = dens) +
  coord_cartesian(xlim = c(-1, 2.5)) +
  xlab("Anomaly") + ylab("Density")
# may1966
# FNN::KL.divergence(total_fits[,809], rnorm(n = 500, mean = 0.4, sd = 0.2))
KLmay1966 = my_KL2(x = total_fits[,809], mu_y = 0.4, sd_y = 0.2)
2*1*KLmay1966
qchisq(c(0.95), df = 2)
hist(rnorm(n = 500, mean = 0.4, sd = 0.2))
abline(v=mean(total_fits[,809]))
# june1966
# FNN::KL.divergence(total_fits[,810], rnorm(n = 500, mean = 1, sd = 0.4))
KLjune1966 = my_KL2(x = total_fits[,810], mu_y = 1, sd_y = 0.4)
2*1*KLjune1966
hist(rnorm(n = 500, mean = 1, sd = 0.4))
abline(v=mean(total_fits[,810]))
# july1966
# FNN::KL.divergence(total_fits[,811], rnorm(n = 500, mean = 0.6, sd = 0.5))
KLjuly1966 = my_KL2(x = total_fits[,811], mu_y = 0.6, sd_y = 0.5)
2*1*KLjuly1966
hist(rnorm(n = 500, mean = 0.6, sd = 0.5))
abline(v=mean(total_fits[,811]))
# august1966
# FNN::KL.divergence(total_fits[,812], rnorm(n = 500, mean = -1.6, sd = 0.4))
KLaugust1966 = my_KL2(x = total_fits[,812], mu_y = -1.6, sd_y = 0.4)
2*1*KLaugust1966
hist(rnorm(n = 500, mean = -1.6, sd = 0.4))
abline(v=mean(total_fits[,812]))

# plot_dat = data.frame(total_fits[,c(789, 809, 810, 811, 812)],
#                       rnorm(n = 2500, mean = 0.9, sd = 0.4),
#                       rnorm(n = 2500, mean = 0.4, sd = 0.2),
#                       rnorm(n = 2500, mean = 1, sd = 0.4),
#                       rnorm(n = 2500, mean = 0.6, sd = 0.5),
#                       rnorm(n = 2500, mean = -1.6, sd = 0.4))
# boxplot(plot_dat[,c(1,6,2,7,3,8,4,9,5,10)], col = 2:3)

comp_data = data.frame(month = c("September 1964", "May 1966", "June 1966", "July 1966", "August 1966"),
                       sat_mean = c(0.9, 0.4, 1, 0.6, -1.6),
                       sat_sd = c(0.4, 0.2, 0.4, 0.5, 0.4),
                       rec_mean = apply(total_fits[,c(789, 809, 810, 811, 812)], 2, mean) %>% round(2),
                       rec_sd = apply(total_fits[,c(789, 809, 810, 811, 812)], 2, sd) %>% round(2),
                       "KL_div*2" = (2 * c(KLseptember1964, KLmay1966, KLjune1966, KLjuly1966, KLaugust1966)) %>% round(2),
                       p_val = 1 - pchisq(2 * c(KLseptember1964, KLmay1966, KLjune1966, KLjuly1966, KLaugust1966), df = 2) %>% round(2))
comp_data
print(xtable::xtable(comp_data), include.rownames=FALSE)

norm_dist = dist(rnorm(n = 100, mean = 0.9, sd = 0.9))
str(norm_dist)
plot(norm_dist)


posterior = function(prior_mean, prior_sd, lik_mean, lik_sd, n) {
  post_mean = prior_sd^2 / (prior_sd^2 + lik_sd^2 / n) * lik_mean +
    lik_sd^2 / (prior_sd^2 + lik_sd^2 / n) * prior_mean
  post_sd = sqrt((1 / prior_sd^2 + n / lik_sd^2)^{-1})
  return(c(post_mean, post_sd))
}
posterior(comp_data[1,2], comp_data[1,3], comp_data[1,4], comp_data[1,5], n = 1) %>% round(2)
posterior(comp_data[2,2], comp_data[2,3], comp_data[2,4], comp_data[2,5], n = 1) %>% round(2)
posterior(comp_data[3,2], comp_data[3,3], comp_data[3,4], comp_data[3,5], n = 1) %>% round(2)
posterior(comp_data[4,2], comp_data[4,3], comp_data[4,4], comp_data[4,5], n = 1) %>% round(2)
posterior(comp_data[5,2], comp_data[5,3], comp_data[5,4], comp_data[5,5], n = 1) %>% round(2)


# ##############################################################################
# analyze the predicted annual cycle over the course of the 20th century
# compute average prediction by sector and month
avg_prediction = data.frame("King_Haakon" = NA, "Ross" = NA,
                            "East_Antarctica" = NA, "Weddell" = NA,
                            "Bellingshausen_Amundsen" = NA)
for (i in 1:dim(fits)[2]) {
  avg_prediction[i, 1:5] = colMeans(fits[,i,])
}
avg_prediction$Total = rowSums(avg_prediction)

# think about these sds
# compute sds
for (i in 1:dim(fits)[2]) {
  avg_prediction[i, 7:11] = apply(fits[,i,], 2, sd)
}
avg_prediction[, 12] = apply(total_fits[], 2, sd)
names(avg_prediction)[7:12] = c("King_Haakon_sd", "Ross_sd", "East_Antarctica_sd", "Weddell_sd",
                                "Bellingshausen_Amundsen_sd", "Total_sd")

# get a prediction for the total
avg_prediction$Month = 1:12
avg_prediction$Season = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1, 1)
avg_prediction$Year = rep(1899:2020, each = 12)
# remove 1899 and 2020
avg_prediction = avg_prediction[!(avg_prediction$Year %in% c(1899, 1980:2020)),]
avg_prediction$Decade = as.character(rep(seq(1900, 1970, by = 10), each = 12*10))

# compute monthly average per decade
decadal_avg = avg_prediction %>%
  group_by(Decade, Month) %>%
  summarize_all(.funs = mean)

decadal_avg$King_Haakon_text = paste0(weights::rd(decadal_avg$King_Haakon, digits=2, add=FALSE), "(",
                                      weights::rd(decadal_avg$King_Haakon_sd, digits=2, add=FALSE), ")")
decadal_avg$Ross_text = paste0(weights::rd(decadal_avg$Ross, digits=2, add=FALSE), "(",
                               weights::rd(decadal_avg$Ross_sd, digits=2, add=FALSE), ")")
decadal_avg$East_Antarctica_text = paste0(weights::rd(decadal_avg$East_Antarctica, digits=2, add=FALSE), "(",
                                          weights::rd(decadal_avg$East_Antarctica_sd, digits=2, add=FALSE), ")")
decadal_avg$Weddell_text = paste0(weights::rd(decadal_avg$Weddell, digits=2, add=FALSE), "(",
                                  weights::rd(decadal_avg$Weddell_sd, digits=2, add=FALSE), ")")
decadal_avg$Bellingshausen_Amundsen_text = paste0(weights::rd(decadal_avg$Bellingshausen_Amundsen, digits=2, add=FALSE), "(",
                                                  weights::rd(decadal_avg$Bellingshausen_Amundsen_sd, digits=2, add=FALSE), ")")
decadal_avg$Total_text = paste0(weights::rd(decadal_avg$Total, digits=2, add=FALSE), "(",
                                weights::rd(decadal_avg$Total_sd, digits=2, add=FALSE), ")")

# King Haakon
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = King_Haakon), show.legend = FALSE) +
  # geom_text(aes(x = Decade, y = Month, label = round(King_Haakon, 2))) +
  geom_text(aes(x = Decade, y = Month, label = King_Haakon_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("King Haakon VII") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_king_haakon.png", device = "png", height = 4, width = 7)
ggplot(decadal_avg) +
  geom_line(aes(group = Decade, col = Decade, x = Month, y = King_Haakon)) +
  viridis::scale_color_viridis("Decade", discrete = TRUE, end = 0.9, direction = -1) +
  ggtitle("King Haakon VII") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_king_haakon_lines.png", device = "png", height = 4, width = 7)


# Ross
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Ross), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Ross_text)) +
  scale_fill_gradient2(limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("Ross Sea") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_ross.png", device = "png", height = 4, width = 7)
ggplot(decadal_avg) +
  geom_line(aes(group = Decade, col = Decade, x = Month, y = Ross)) +
  viridis::scale_color_viridis("Decade", discrete = TRUE, end = 0.9, direction = -1) +
  ggtitle("Ross Sea") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_ross_lines.png", device = "png", height = 4, width = 7)

# East Antarctica
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = East_Antarctica), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = East_Antarctica_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  ggtitle("East Antarctica") +
  scale_y_continuous(breaks = 1:12) +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_east_antarctica.png", device = "png", height = 4, width = 7)
ggplot(decadal_avg) +
  geom_line(aes(group = Decade, col = Decade, x = Month, y = East_Antarctica)) +
  viridis::scale_color_viridis("Decade", discrete = TRUE, end = 0.9, direction = -1) +
  ggtitle("East Antarctica") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_east_antarctica_lines.png", device = "png", height = 4, width = 7)


# Weddell Sea
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Weddell), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Weddell_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.4, 0.4), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  ggtitle("Weddell Sea") +
  scale_y_continuous(breaks = 1:12) +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_weddell.png", device = "png", height = 4, width = 7)
ggplot(decadal_avg) +
  geom_line(aes(group = Decade, col = Decade, x = Month, y = Weddell)) +
  viridis::scale_color_viridis("Decade", discrete = TRUE, end = 0.9, direction = -1) +
  ggtitle("Weddell Sea") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_weddell_lines.png", device = "png", height = 4, width = 7)

# Bellingshausen Amundsen
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Bellingshausen_Amundsen), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Bellingshausen_Amundsen_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.2, 0.2), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("Bellingshausen Amundsen Sea") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_bellingshausen_amundsen.png", device = "png", height = 4, width = 7)
ggplot(decadal_avg) +
  geom_line(aes(group = Decade, col = Decade, x = Month, y = Bellingshausen_Amundsen)) +
  viridis::scale_color_viridis("Decade", discrete = TRUE, end = 0.9, direction = -1) +
  ggtitle("Bellingshausen Amundsen Sea") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_bellingshausen_amundsen_lines.png", device = "png", height = 4, width = 7)


# Total
ggplot(decadal_avg) +
  geom_raster(aes(x = Decade, y = Month, fill = Total), show.legend = FALSE) +
  geom_text(aes(x = Decade, y = Month, label = Total_text)) +
  scale_fill_gradient2("Anomaly",
                       limits = c(-0.4, 0.4), oob = scales::squish,
                       low = "firebrick",
                       mid = "white",
                       high = "steelblue") +
  scale_y_continuous(breaks = 1:12) +
  ggtitle("Total") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_total.png", device = "png", height = 4, width = 7)
ggplot(decadal_avg) +
  geom_line(aes(group = Decade, col = Decade, x = Month, y = Total)) +
  viridis::scale_color_viridis("Decade", discrete = TRUE, end = 0.9, direction = -1) +
  ggtitle("Total") +
  theme_minimal()
ggsave("../plots/BayesianModel/decadal_avg_total_lines.png", device = "png", height = 4, width = 7)

# # summarize by season
# decadal_avg_s = avg_prediction %>%
#   group_by(Decade, Season) %>%
#   summarize_all(.funs = mean)
#
# # Total
# ggplot(decadal_avg_s) +
#   geom_raster(aes(x = Decade, y = Season, fill = Total)) +
#   scale_fill_gradient2(low = "firebrick",
#                        mid = "white",
#                        high = "steelblue") +
#   scale_y_continuous(breaks = 1:4) +
#   theme_minimal()






