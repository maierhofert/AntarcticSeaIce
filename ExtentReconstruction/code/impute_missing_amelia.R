library("Amelia")
source("read_Fogt_predictor_data.R")
# amelia does not impute entirely missing variables
# that's ok, we have at least one covariate for every day


# fit imputation model
data_ghcnd$ts = 1:nrow(data_ghcnd)
data_ghcnd$cs = 1

# # get longer leads and lags
# datASN2011$sTMIN_lag3 = lag(datASN2011$sTMIN, n = 3, default = NA)
# datASN2011$sTMIN_lag7 = lag(datASN2011$sTMIN, n = 7, default = NA)
# datASN2011$sTMIN_lead3 = lead(datASN2011$sTMIN, n = 3, default = NA)
# datASN2011$sTMIN_lead7 = lead(datASN2011$sTMIN, n = 7, default = NA)

measure_vars = ! (names(data_ghcnd) %in% c("Year", "Month", "Day", "Date", "Doy"))
amelia_mod = amelia(x = data_ghcnd[, measure_vars],
                    ts = "ts", polytime = NULL,
                    cs = "cs",
                    parallel = "multicore", ncpus = 10,
                    # lags = c("sTMIN", "sTMAX"),
                    # leads = c("sTMIN", "sTMAX"),
                    m = 1)
summary(amelia_mod)
summary(amelia_mod$imputations)
amelia_mod$covMatrices %>% round(2)
# plot(amelia_mod)

# plot imputations
tscsPlot(amelia_mod, cs = 1, var = "sTMIN")
         
# plot the actual observations
plot(datASN2011$sTMIN[450:550], pch = 16)
lines(amelia_mod$imputations$imp1$sTMIN, type = "l", col = 2)
lines(amelia_mod$imputations$imp2$sTMIN, type = "l", col = 3)
lines(amelia_mod$imputations$imp3$sTMIN, type = "l", col = 4)

plot(amelia_mod$imputations$imp1$sTMAX, type = "l")
lines(amelia_mod$imputations$imp2$sTMAX, type = "l", col = 2)

# Some observations are missing in both variables
# this should not be a problem when using all covariates
table(is.na(datASN2011$TMIN[22000:23000]), is.na(datASN2011$TMAX[22000:23000]))
table(is.na(datASN2011$TMIN), is.na(datASN2011$TMAX))

# create a full model
# still really fast
amelia_mod_full = amelia(x = datASN2011[, c("ts", "cs", "sTMIN", "sTMAX")],
                    ts = "ts", polytime = NULL,
                    cs = "cs",
                    lags = c("sTMIN", "sTMAX"),
                    leads = c("sTMIN", "sTMAX"),
                    m = 5)
plot(amelia_mod_full$imputations$imp1$sTMIN[450:550], type = "l")

plot(amelia_mod_full$imputations$imp1$sTMIN[800:900], type = "l")
lines(amelia_mod_full$imputations$imp1$sTMIN[800:900], type = "p", pch = 16, 
     col = is.na(datASN2011$TMIN[800:900]) + 2)
