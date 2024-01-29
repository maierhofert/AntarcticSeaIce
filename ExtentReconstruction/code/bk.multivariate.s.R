# multivariate version of bk.1.s

library(bayesbackcast)
library(rstan)
fname <- c("01", "02", "03", "04", "05")
n_cores <- 4
chains <- 4
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
all_dat = readRDS("../data/complete_dataset/all_data.RDS")
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
#odd <- rep(c(TRUE, FALSE), length.out=length(buff))
#buff <- all_dat$tdate >= all_dat$tdate[buff.start] & all_dat$tdate < 2021.5 & odd
all_dat <- all_dat[buff,]
all_x <- all_x[buff,]
all_dat$sie <- all_dat[,sel_region+3]
# all sea ice is observed
sum(is.na(all_dat$sie))
all_dat$sie[is.na(all_dat$sie)] <- -1000
sie <- ts(data = all_dat$sie, start = min(all_dat$tdate), end = max(all_dat$tdate),
          frequency = 12)

# Add year indicators (not currently used)
#year_x <- model.matrix(~ -1 + factor(all_dat$year))[,-1]
#year_x <- year_x[,ncol(year_x)-c(5,0),drop=FALSE]

# get data
#start.date <- which.max(all_dat$tdate > 1978.82)
#start.date
start <- all_dat$tdate > 1979.05
all_dat$tdate[start][1:10]
no1978 <- seq_along(sie[,1])[start]
table(all_dat$year[no1978])
all_x = all_x[no1978,]
sie[all_dat$tdate > 1987.91 & all_dat$tdate < 1988] <- -1000
sie[all_dat$tdate > 1987.55 & all_dat$tdate < 1987.7] <- -1000
sum(sie[no1978] < -999)

# # Add splines
# p <- ncol(all_x)
# p
# call_x = colnames(all_x)
#all_x = cbind(all_x,all_x,all_x,all_x)
#Doy <- rep(1:365,nrow(all_x)/365)
#for(j in 2:4){
#  colnames(all_x)[(1:p)+(j-1)*p] <- paste0("b",j-1,"_",call_x)
#  all_x[,(1:p)+(j-1)*p] <- sweep(all_x[,(1:p)+(j-1)*p],1,pmsMat[Doy,j-1],"*")
#}

p <- ncol(all_x)
n <- nrow(all_x)
p_nonzero <- 8 # prior guess for the number of relevant variables
p_nonzero <- 15

# select number of retained pc here
all_pc <- prcomp(all_x, center=TRUE, scale=TRUE)$x[,1:100]

# number of autregressive components
Kar <- 4
# 0 moving average components
Kma <- 0
# code up an undifferenced version
mdl <- arima(ts = sie, order = c(Kar,0,Kma), xreg = all_pc,
             slab_scale=0.025, tprior=FALSE,
             series.name = "all",seldat=no1978,
             lagAR=c(1,2,3,12),lagMA=0, p_nonzero=p_nonzero)
str(mdl)

# library(roxygen2)
# pkgbuild::compile_dll(path = ".", force = T)
# # should check and only recompile new files
# pkgbuild::compile_dll(path = ".", force = F)
# or
# pkgbuild::compile_dll(path = ".")
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
  fit <- varstan(model = mdl, iter=30, chains=1, tree.depth=10)
  fit <- varstan(mode = mdl, iter=1000, chains=chains, tree.depth=10)
  # saveRDS(fit, file = paste0("../output/bc.",fname,".od.",vname,".rds"))
  # save location when in bayesbackcast package
  saveRDS(fit, file = paste0("../RECONSTRUCTION/output/bc.",fname,".od.",vname, ".nz.", p_nonzero, ".rds"))
}else{
  fit <- readRDS(file = paste0("../output/bc.",fname,".od.",vname,".rds"))
}
# pdf(paste0("../output/f.",fname,".od.",vname,".pdf"))
# check_residuals(fit)
# #autoplot(forecast(object = fit,h = 12))
# dev.off()

# str(fit)
library("dplyr")
fit$stanfit@sim %>% str()
extract(fit$stanfit) %>% str()
# look at parameter summary
bayesbackcast:::summary.varstan(fit, pars = c("beta0C", "brg"))
# reformat to 5x5
bayesbackcast:::summary.varstan(fit, pars = "Sigma")

bayesbackcast:::summary.varstan(fit, pars = "ar")
eps = bayesbackcast:::summary.varstan(fit, pars = "epsilon")
head(eps)
tail(eps)
eps %>%
  mutate(sector = rep(1:5, each = 498)) %>%
  group_by(sector) %>%
  summarise_all(mean, na.rm = TRUE)

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

pars1 <- c("ar[1,1]","ma[1,1]","ar[1,2]","ma[1,2]","sigma0","beta_0")
pars1 <- c("ar[1,1]","ma[1,1]","sigma0","beta_0")
pars1 <- c("ar[1,1]","sigma0","beta_0")
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

if(F){
  loo1 <- loo(fit, cores = 8)
  pdf(paste0("../output/g.",fname,".od.",vname,".pdf"))
  plot(loo1)
  print(loo1)
}

if(any(is.na(all_dat$sie))){
  ypost = (rstan::extract(fit$stanfit, "latent_data", permuted = FALSE))
  y <- NULL
  for(i in dim(ypost)[2]){
    y <- rbind(y, ypost[,i,])
  }
  yactual <- fit$model$yall
  miss <- yactual < -999
  #missp <- miss & c(FALSE,miss[-length(miss)]) & c(miss[-1],FALSE)
  yactual[miss] <- NA
  #yactual <- yactual[-1]
  sum(miss)
  ylatent <- yactual
  ylatent[miss] <- y[1,]
  plot(x=all_dat$tdate[miss],y=ylatent[miss],ylim=c(-1.5,1.5),xlim=c(1987,1989.4),xlab="Date",
       ylab="anomoly", col=3,pch=16, type='n', main="Stochastic reconstruction of satellite era",
       sub="green are stochastically imputed")
  #for(i in 1:dim(ypost)[1]){
  for(i in dim(ypost)[1]){
    ylatent <- yactual
    ylatent[miss] <- y[i,]
    ylatent[!miss] <- NA
    ylatent <- ylatent[-1]
    #points(x=all_dat$tdate[miss],y=ylatent[miss],xlim=c(1987,1989.4),col=3,pch=16,cex=0.6)
    lines(x=all_dat$tdate,y=yactual[-1],xlim=c(1987,1989.4),col=1)
    lines(x=all_dat$tdate,y=ylatent,xlim=c(1987,1989.4),col=3)
  }
  plot(  x=all_dat$tdate[ miss],y=ylatent[miss],ylim=c(-1.5,1.5),xlim=c(1987,1989.4),xlab="Date",
         ylab="anomoly", col=3,pch=16, cex=0.2, main="Stochastic reconstruction of satellite era")
  lines(x=all_dat$tdate,y=ylatent,xlim=c(1987,1989.4),col=3)
}

dev.off()
