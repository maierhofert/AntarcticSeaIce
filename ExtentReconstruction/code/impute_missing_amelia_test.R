library("Amelia")
# amelia does not impute entirely missing variables
# that should be ok, we should have at least one covariate for every day


# fit imputation model
datASN2011$ts = 1:nrow(datASN2011)
datASN2011$cs = 1

# amelia_mod = amelia(x = datASN2011[450:550, c("ts", "cs", "TMIN", "TMAX")],
#                     ts = "ts", polytime = 3,
#                     cs = "cs",
#                # x = datASN2011[, c("TMIN", "TMAX")], 
#                # x = datASN2011[22000:23000, c("TMIN", "TMAX")],
#               m = 5)

# get longer leads and lags
datASN2011$sTMIN_lag3 = lag(datASN2011$sTMIN, n = 3, default = NA)
datASN2011$sTMIN_lag7 = lag(datASN2011$sTMIN, n = 7, default = NA)
datASN2011$sTMIN_lead3 = lead(datASN2011$sTMIN, n = 3, default = NA)
datASN2011$sTMIN_lead7 = lead(datASN2011$sTMIN, n = 7, default = NA)


amelia_mod = amelia(x = datASN2011[450:550, c("ts", "cs", "sTMIN", "sTMAX", 
                                              "sTMIN_lag3", "sTMIN_lag7",
                                              "sTMIN_lead3", "sTMIN_lead7")],
                    ts = "ts", polytime = NULL,
                    cs = "cs",
                    # lags = c("sTMIN", "sTMAX"),
                    # leads = c("sTMIN", "sTMAX"),
                    m = 5)
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
