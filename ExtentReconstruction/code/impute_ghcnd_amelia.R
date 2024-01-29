library("Amelia")
source("read_Fogt_predictor_data_ghcnd.R")
# amelia does not impute entirely missing variables
# that's ok, we have at least one covariate for every day

# subset to data after 1900
data_ghcnd_detrended = data_ghcnd_detrended[data_ghcnd_detrended$Date >= as.Date("1900-01-01"),]
# subset to data where first observed observation is before 1979
first_observed = apply(data_ghcnd_detrended, 2, function (x) {
  min(which(!is.na(x))[1], length(x), na.rm = TRUE)
})
first_observed_date = data_ghcnd_detrended$Date[first_observed]
delete = first_observed_date > as.Date("1979-01-01")
# how many columns are affected
sum(delete)
# do deletion
data_ghcnd_detrended = data_ghcnd_detrended[, !delete]

# add linear time to data
data_ghcnd_detrended$ts = 1:nrow(data_ghcnd_detrended)
data_ghcnd_detrended$cs = 1
# data_ghcnd_detrended$white_noise = rnorm(nrow(data_ghcnd_detrended))

# split into groups of similar covariates
idvars = c("ts", "cs", "Year", "Month", "Day", "Date", "Doy")
group = substr(names(data_ghcnd_detrended), 6, 8)
group[group == "AYW"] = "AYM"
group[names(data_ghcnd_detrended) %in% idvars] = "idvar"
# group[names(data_ghcnd_detrended) == "white_noise"] = "white_noise"
table(group)

# split the ASN group into 4
group[4:21] = "ASN1"
group[22:55] = "ASN2"
group[56:73] = "ASN3"
group[74:92] = "ASN4"
table(group)

# # check group membership
# names(data_ghcnd_detrended)[group == "ASN1"]
# names(data_ghcnd_detrended)[group == "ASN2"]
# names(data_ghcnd_detrended)[group == "ASN3"]
# names(data_ghcnd_detrended)[group == "ASN4"]

# do imputation
# ##############################################################################
# for group of ids in ASN1
namesASN1 = names(data_ghcnd_detrended)[group == "ASN1"]

Sys.time()
amelia_mod_ASN1 = amelia(x = data_ghcnd_detrended[, group %in% c("ASN1", "idvar")],
                    idvars = idvars[-(1:2)], # take out ts and cs from idvars
                    ts = "ts", polytime = NULL, 
                    cs = "cs",
                    parallel = "multicore", ncpus = 10,
                    # lags = names6,
                    # leads = names6,
                    m = 2)
Sys.time()
summary(amelia_mod_ASN1)
summary(amelia_mod_ASN1$imputations)

# # plot imputations
# for (i in 1:length(namesASN1)) {
#   tscsPlot(amelia_mod_ASN1, cs = 1, var = namesASN1[i], main = namesASN1[i])
# }
saveRDS(amelia_mod_ASN1, "../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN1.RDS")

# ##############################################################################
# for group of ids in ASN2
namesASN2 = names(data_ghcnd_detrended)[group == "ASN2"]
Sys.time()
amelia_mod_ASN2 = amelia(x = data_ghcnd_detrended[, group %in% c("ASN2", "idvar")],
                         idvars = idvars[-(1:2)], # take out ts and cs from idvars
                         ts = "ts", polytime = NULL, 
                         cs = "cs",
                         parallel = "multicore", ncpus = 10,
                         # lags = names6,
                         # leads = names6,
                         m = 2)
Sys.time()
summary(amelia_mod_ASN2)
summary(amelia_mod_ASN2$imputations)

# # plot imputations
# for (i in 1:length(namesASN2)) {
#   tscsPlot(amelia_mod_ASN2, cs = 1, var = namesASN2[i], main = namesASN2[i])
# }
saveRDS(amelia_mod_ASN2, "../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN2.RDS")


# ##############################################################################
# for group of ids in ASN3
namesASN3 = names(data_ghcnd_detrended)[group == "ASN3"]
Sys.time()
amelia_mod_ASN3 = amelia(x = data_ghcnd_detrended[, group %in% c("ASN3", "idvar")],
                         idvars = idvars[-(1:2)], # take out ts and cs from idvars
                         ts = "ts", polytime = NULL, 
                         cs = "cs",
                         parallel = "multicore", ncpus = 10,
                         # lags = names6,
                         # leads = names6,
                         m = 2)
Sys.time()
summary(amelia_mod_ASN3)
summary(amelia_mod_ASN3$imputations)

# # plot imputations
# for (i in 1:length(namesASN3)) {
#   tscsPlot(amelia_mod_ASN3, cs = 1, var = namesASN3[i], main = namesASN3[i])
# }
saveRDS(amelia_mod_ASN3, "../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN3.RDS")

# ##############################################################################
# for group of ids in ASN4
namesASN4 = names(data_ghcnd_detrended)[group == "ASN4"]
Sys.time()
amelia_mod_ASN4 = amelia(x = data_ghcnd_detrended[, group %in% c("ASN4", "idvar")],
                         idvars = idvars[-(1:2)], # take out ts and cs from idvars
                         ts = "ts", polytime = NULL, 
                         cs = "cs",
                         parallel = "multicore", ncpus = 10,
                         # lags = names6,
                         # leads = names6,
                         m = 2)
Sys.time()
summary(amelia_mod_ASN4)
summary(amelia_mod_ASN4$imputations)

# # plot imputations
# for (i in 1:length(namesASN4)) {
#   tscsPlot(amelia_mod_ASN4, cs = 1, var = namesASN4[i], main = namesASN4[i])
# }
saveRDS(amelia_mod_ASN4, "../data/fogt_predictor_data_imputed/amelia_ghcnd_groupASN4.RDS")

# ##############################################################################
# for group of ids in AYM and AYW
namesAYM = names(data_ghcnd_detrended)[group == "AYM"]
Sys.time()
amelia_mod_AYM = amelia(x = data_ghcnd_detrended[, group %in% c("AYM", "idvar", )],
                         idvars = idvars[-(1:2)], # take out ts and cs from idvars
                         ts = "ts", polytime = NULL, 
                         cs = "cs",
                         parallel = "multicore", ncpus = 10,
                         # lags = names6,
                         # leads = names6,
                         m = 2)
Sys.time()
summary(amelia_mod_AYM)
summary(amelia_mod_AYM$imputations)

# barely goes back to 1950s
image(x = data_ghcnd$Date, y = 1:length(namesAYM), 
      z = miss_meas[,colnames(miss_meas) %in% namesAYM], 
      xlab = "Day", ylab = "Station", 
      col = c("black", "white"))

# plot imputations
for (i in 1:length(namesAYM)) {
  tscsPlot(amelia_mod_AYM, cs = 1, var = namesAYM[i], main = namesAYM[i])
}
saveRDS(amelia_mod_AYM, "../data/fogt_predictor_data_imputed/amelia_ghcnd_groupAYM.RDS")


