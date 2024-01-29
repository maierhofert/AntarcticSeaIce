library("Amelia")
source("read_Fogt_predictor_data_monthly.R")
# amelia does not impute entirely missing variables
# that's ok, we have at least one covariate for every day

# subset to data after 1899
data_mon_detrended = data_mon_detrended[between(data_mon_detrended$Date, as.Date("1899-01-01"), as.Date("2021-01-01")),]
# subset to data where first observed observation is before 1979
first_observed = apply(data_mon_detrended, 2, function (x) {
  min(which(!is.na(x))[1], length(x), na.rm = TRUE)
})
first_observed_date = data_mon_detrended$Date[first_observed]

# is there data left that starts after the beginning of the satellite era
delete = first_observed_date > as.Date("1979-01-01")
# how many columns are affected
sum(delete)
# do deletion
data_mon_detrended = data_mon_detrended[, !delete]

# add linear time to data
data_mon_detrended$ts = 1:nrow(data_mon_detrended)
data_mon_detrended$cs = 1
# data_mon_detrended$white_noise = rnorm(nrow(data_mon_detrended))

# split into groups of similar covariates
idvars = c("ts", "cs", "Year", "Month", "Date", "Doy")
group = substr(names(data_mon_detrended), 1, 1)
group[names(data_mon_detrended) %in% idvars] = "idvar"
# group[names(data_mon_detrended) == "white_noise"] = "white_noise"

# split the group into two
group2 = substr(names(data_mon_detrended), 1, 2)
group[group2 %in% c("83", "85", "87")] = "8" # 83
group[group2 %in% c("88", "89")] = "8" # 88
group[group2 %in% c("Ma", "Sc")] = "8" # 88
group[group2 %in% c("IP", "PD", "SA", "AM", "SO", "SS")] = 9 # "ClimateInd"

# how many variables per group
table(group)

# apply(data_mon_detrended, 2, function(x) {
#   sum(is.na(x))
# })

# do imputation

# for group of ids starting with 6
names6 = names(data_mon_detrended)[group == "6"]
# less than 2 minutes
Sys.time()
amelia_mod_6 = amelia(x = data_mon_detrended[, group %in% c("6", "idvar")],
                      idvars = idvars[-(1:2)], # take out ts and cs from idvars
                      ts = "ts", polytime = NULL, 
                      cs = "cs",
                      parallel = "multicore", ncpus = 10,
                      lags = names6,
                      # leads = names6,
                      m = 5)
Sys.time()
summary(amelia_mod_6)
summary(amelia_mod_6$imputations)

# # plot imputations
# for (i in 1:length(names6)) {
#   pdf(paste0("../plots/covariates/imputations/", names6[i], ".pdf"), height = 7, width = 10)
#   tscsPlot(amelia_mod_6, cs = 1, var = names6[i], main = names6[i])
#   dev.off()
# }
saveRDS(amelia_mod_6, "../data/fogt_predictor_data_imputed/amelia_monthly_group6.RDS")

# for group of ids starting with 8
names8 = names(data_mon_detrended)[group == "8"]
Sys.time()
amelia_mod_8 = amelia(x = data_mon_detrended[, group %in% c("8", "idvar")],
                       idvars = idvars[-(1:2)], # take out ts and cs from idvars
                       ts = "ts", polytime = NULL, 
                       cs = "cs",
                       parallel = "multicore", ncpus = 10,
                       lags = names8,
                       # leads = names83,
                       m = 5)
Sys.time()
summary(amelia_mod_8)
summary(amelia_mod_8$imputations)
# # plot imputations
# for (i in 1:length(names8)) {
#   pdf(paste0("../plots/covariates/imputations/", names8[i], ".pdf"), height = 7, width = 10)
#   tscsPlot(amelia_mod_8, cs = 1, var = names8[i], main = names8[i])
#   dev.off()
# }
saveRDS(amelia_mod_8, "../data/fogt_predictor_data_imputed/amelia_monthly_group8.RDS")



# # for group of ids starting with 83
# names83 = names(data_mon_detrended)[group == "83"]
# Sys.time()
# amelia_mod_83 = amelia(x = data_mon_detrended[, group %in% c("83", "idvar")],
#                        idvars = idvars[-(1:2)], # take out ts and cs from idvars
#                        ts = "ts", polytime = NULL, 
#                        cs = "cs",
#                        parallel = "multicore", ncpus = 10,
#                        lags = names83,
#                        # leads = names83,
#                        m = 5)
# Sys.time()
# summary(amelia_mod_83)
# summary(amelia_mod_83$imputations)
# # # plot imputations
# # for (i in 1:length(names83)) {
# #   tscsPlot(amelia_mod_83, cs = 1, var = names83[i], main = names83[i])
# # }
# saveRDS(amelia_mod_83, "../data/fogt_predictor_data_imputed/amelia_monthly_group83.RDS")
# 
# # # for group of ids starting with 88
# # # the problem here is that there is no data available for the first 39 months,
# # # which is about 3 years
# # data88 = data_mon_detrended[, group %in% c("88")]
# # all_miss = apply(data88, 1, function(x) {
# #   mean(is.na(x)) == 1
# # })
# # table(all_miss)
# # plot(all_miss)
# 
# names88 = names(data_mon_detrended)[group == "88"]
# Sys.time()
# amelia_mod_88 = amelia(x = data_mon_detrended[, group %in% c("83", "88", "idvar", "white_noise")],
#                        idvars = idvars[-(1:2)], # take out ts and cs from idvars
#                        ts = "ts", polytime = NULL, 
#                        cs = "cs",
#                        parallel = "multicore", ncpus = 10,
#                        lags = names88,
#                        # leads = names88,
#                        m = 5)
# Sys.time()
# summary(amelia_mod_88)
# summary(amelia_mod_88$imputations)
# # plot imputations
# for (i in 1:length(names88)) {
#   tscsPlot(amelia_mod_88, cs = 1, var = names88[i], main = names88[i])
# }
# saveRDS(amelia_mod_88, "../data/fogt_predictor_data_imputed/amelia_monthly_group88.RDS")

# for group of ids starting with 9
Sys.time()
names9 = names(data_mon_detrended)[group == "9"]
amelia_mod_9 = amelia(x = data_mon_detrended[, group %in% c("9", "idvar")],
                      idvars = idvars[-(1:2)], # take out ts and cs from idvars
                      ts = "ts", polytime = NULL, 
                      cs = "cs",
                      parallel = "multicore", ncpus = 10,
                      lags = names9,
                      # leads = names9,
                      m = 5)
Sys.time()
summary(amelia_mod_9)
summary(amelia_mod_9$imputations)
# plot imputations
for (i in 1:length(names9)) {
  pdf(paste0("../plots/covariates/imputations/", names9[i], ".pdf"), height = 7, width = 10)
  tscsPlot(amelia_mod_9, cs = 1, var = names9[i], main = names9[i])
  dev.off()
}
saveRDS(amelia_mod_9, "../data/fogt_predictor_data_imputed/amelia_monthly_group9.RDS")

