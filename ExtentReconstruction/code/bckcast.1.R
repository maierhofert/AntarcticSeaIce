library(bayesbackcast)
library(rstan)
# fname <- "01"
# vname <- "s"
# n_cores <- 4
# chains <- 4
# sel_region <- match(fname,c("01","02","03","04","05"))
# fname <- "01.all"
# name_region <- c("King Hakon VII", "Ross", "East Antarctica",
#                  "Weddell", "Bellingshausen-Amundsen Sea")[sel_region]
# rstan_options(auto_write = TRUE, cores = n_cores)
# options(mc.cores=n_cores)
# library(bayesplot)
# library(loo)
# set.seed(1)
#
#
# # read in data
# all_dat = readRDS("../data/complete_dataset/all_data.RDS")
# # subset to before 2021
# all_dat = all_dat[all_dat$Year <= 2020,]
# # # subset to first 100 PC
# # all_dat = all_dat[,1:109]
# all_dat$tdate = all_dat$Year + all_dat$Month / 12 - 0.5
# # extASN <- grep("ASN0",colnames(all_dat),fixed=T)
# # all_dat <- all_dat[,-extASN]
# #inc <- colnames(all_dat)[c(1:8,grep("ppb",colnames(all_dat),fixed=TRUE))]
# #inc <- c(inc, "VMarsh_O_Higgins_SLP","VMarsh_O_Higgins_tmp","VScott_Base_mcMurdo_SLP","VMarsh_O_Higgins_tmp")
# #inc <- inc[-grep("_l7",inc,fixed=T)]
# #inc
# #all_dat <- all_dat[,inc]
# all_x <- all_dat[,-c(1:9)]
# #round(cor(all_dat[,inc],use="complete"),2)
# #buff.start <- which.max(all_dat$tdate > 1978.82)-1
# #buff <- all_dat$tdate >= all_dat$tdate[buff.start] & all_dat$tdate < 2021.5
# buff <- all_dat$Date >= as.Date("1979-01-01") & all_dat$Date < as.Date("2021-07-30")
# #odd <- rep(c(TRUE, FALSE), length.out=length(buff))
# #buff <- all_dat$tdate >= all_dat$tdate[buff.start] & all_dat$tdate < 2021.5 & odd
# all_dat <- all_dat[buff,]
# all_x <- all_x[buff,]
# all_dat$sie <- all_dat[,sel_region+3]
# sum(is.na(all_dat$sie))
# all_dat$sie[is.na(all_dat$sie)] <- -1000
# sie <- ts(data = all_dat$sie, start = min(all_dat$tdate), end = max(all_dat$tdate),
#           frequency = 12)
#
# # Add year indicators (not currently used)
# #year_x <- model.matrix(~ -1 + factor(all_dat$year))[,-1]
# #year_x <- year_x[,ncol(year_x)-c(5,0),drop=FALSE]
#
# # Add splines
# p <- ncol(all_x)
# p
# call_x = colnames(all_x)
# #all_x = cbind(all_x,all_x,all_x,all_x)
# #Doy <- rep(1:365,nrow(all_x)/365)
# #for(j in 2:4){
# #  colnames(all_x)[(1:p)+(j-1)*p] <- paste0("b",j-1,"_",call_x)
# #  all_x[,(1:p)+(j-1)*p] <- sweep(all_x[,(1:p)+(j-1)*p],1,pmsMat[Doy,j-1],"*")
# #}
#
# p <- ncol(all_x)
# n <- nrow(all_x)
# p_nonzero <- 8 # prior guess for the number of relevant variables
#
# fit <- readRDS(file = paste0("bc.",fname,".od.",vname,".rds"))
# paste0("bc.",fname,".od.",vname,".rds")
#
# all_pc <- prcomp(all_x, center=TRUE, scale=TRUE)$x[,1:100]

#debug(backcast)
dim(all_x)
dim(fit$model$X)
pp = posterior_predict(object = fit, h = nrow(all_pc), newdata = all_pc, robust = FALSE,
                       draws = 10, seed = NULL, backcast=TRUE,
                       newtdate=all_dat$tdate, no.volatility=FALSE)
#saveRDS(fc, file = paste0("bkcast.",fname,".rds"))
#
dim(all_dat)
dim(pp)
pdf(paste0("bkcast.",fname,".pdf"))
plot(x=all_dat$tdate, y=pp[1,], xlab="date", ylab="anomoly",
     main=paste0("Reconstructed anomoly for ",name_region), type='l')
length(all_dat$tdate)
plot(x=all_dat$tdate, y=pp[1,], xlab="date", ylab="anomoly", xlim=c(1978,1980),
     main=paste0("Reconstructed anomoly for ",name_region), type='l')
fc <- backcast(fit, newdata=all_pc, newtdate=all_dat$tdate, no.volatility=FALSE)
autoplot(fc, main=paste0("Reconstruction and data for ",name_region))



