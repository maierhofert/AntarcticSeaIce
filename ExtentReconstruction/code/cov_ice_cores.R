# # read ice core data collection from Thomas et al (2019)
# # this has just an annual resolution
# # this is just snow accumulation and no chemicals
# ice_core_dat = readxl::read_xlsx("../data/ice_cores/Ant2k_RegionalComposites_Thomas_2017_Dec_v2.xlsx",
#                                  sheet = "Original data", skip = 6)
# 
# head(ice_core_dat)
# ice_core_dat_sources = readxl::read_xlsx("../data/ice_cores/Ant2k_RegionalComposites_Thomas_2017_Dec_v2.xlsx",
#                                          sheet = "Data Sources", skip = 4)
# head(ice_core_dat_sources)
# ice_core_dat_regional = readxl::read_xlsx("../data/ice_cores/Ant2k_RegionalComposites_Thomas_2017_Dec_v2.xlsx",
#                                           sheet = "Regional Composites", skip = 7)
# head(ice_core_dat_regional)



# read ice core data
# this one has a sub-annual resolution which is awesome
SPC14 = read.csv("../data/ice_cores/SPC14_Ion_Chemistry_2020.csv")
head(SPC14)

# filter our rows where all observations are missing
SPC14 = filter(SPC14, !(is.na(Chloride..ppb.) & is.na(Nitrate..ppb.) & is.na(Sulfate..ppb.) &
                          is.na(Sodium..ppb.) & is.na(Magnesium..ppb.) & is.na(Calcium..ppb.)))

SPC14$avgyear_before1950 = rowMeans(SPC14[,c("Bottom.Age..Years.before.1950.", "Top.Age..years.before.1950.")])
SPC14$Doy = (SPC14$avgyear_before1950 %% 1) * 365.2422

SPC14$DayBeforeJan1_1950 = (SPC14$avgyear_before1950 * 365.2422) %>% round(0)
SPC14$Date = as.Date(-SPC14$DayBeforeJan1_1950, origin = as.Date("1950-01-01"))
summary(SPC14$Date)

# filter to after 1900
SPC14 = SPC14[SPC14$Date >= as.Date("1899-01-01"),]
summary(SPC14)

hist(as.numeric(diff.Date(SPC14$Date)))
plot(SPC14$Date, c(0, as.numeric(diff.Date(SPC14$Date))))


# detrend data
SPC14 = SPC14 %>%
  mutate(Doy = lubridate::yday(Date))

SPC14_detrended = SPC14
measure_vars = names(SPC14)[5:10]
for (i in 1:length(measure_vars)) {
  this.var = measure_vars[i]
  this.gam = mgcv::gam(formula(paste0(this.var, " ~ s(Doy, bs = 'cc')")), data = SPC14)
  # anova(this.gam)
  # plot(this.gam)
  SPC14_detrended[,this.var] = SPC14[,this.var] - unlist(predict(this.gam, newdata = SPC14))
  SPC14_detrended[,this.var] = c(scale(SPC14_detrended[,this.var], center = TRUE, scale = TRUE))
}
summary(SPC14_detrended)

# impute missing values
library("Amelia")
# add linear time to data
SPC14_detrended$ts = 1:nrow(SPC14_detrended)
SPC14_detrended$cs = 1
# SPC14_detrended$white_noise = rnorm(nrow(SPC14_detrended))

idvars = names(SPC14_detrended)[!(names(SPC14_detrended) %in% measure_vars)]

# instantaneous
Sys.time()
amelia_mod_iceCores = amelia(x = SPC14_detrended,
                      idvars = idvars[-(9:10)], # take out ts and cs from idvars
                      ts = "ts", polytime = NULL, 
                      cs = "cs",
                      parallel = "multicore", ncpus = 10,
                      lags = measure_vars,
                      leads = measure_vars,
                      m = 5)
Sys.time()
summary(amelia_mod_iceCores)
tscsPlot(amelia_mod_iceCores, "Calcium..ppb.", cs = 1)


# do a smooth interpolation
# smooth interpolation
Dates <- seq.Date(lubridate::ymd("1899-01-01"), lubridate::ymd("2014-12-31"), by = 1)
SPC14_daily = data.frame(Date = Dates)



# TODO fix this in monthly data as well
for (i in 5:10) {
  smoothxy = ksmooth(amelia_mod_iceCores[[1]][["imp1"]]$Date, 
                     amelia_mod_iceCores[[1]][["imp1"]][,i], 
                     kernel = "normal", bandwidth = 30.4)
  # plot(SPC14$Date, SPC14[,i], type = "l")
  # plot(smoothxy$x, smoothxy$y, type = "l")
  SplineFun <- splinefun(x = smoothxy$x, y = smoothxy$y, method = "natural")
  SplineFit <- SplineFun(Dates)
  plot(Dates, SplineFit, type = "l")
  SPC14_daily[, names(SPC14)[i]] = SplineFit
}
head(SPC14_daily)
iceCore_cov_daily = SPC14_daily

# # plot the data
# plot(SPC14$Date, SPC14$Chloride..ppb., type = "l")
# plot(SPC14_daily$Date, SPC14_daily$Chloride..ppb., type = "l")
# 
# plot(SPC14$Date, SPC14$Nitrate..ppb., type = "l")
# plot(SPC14_daily$Date, SPC14_daily$Nitrate..ppb., type = "l")
# 
# plot(SPC14$Date, SPC14$Sulfate..ppb., type = "l")
# plot(SPC14_daily$Date, SPC14_daily$Sulfate..ppb., type = "l")
# 
# plot(SPC14$Date, SPC14$Sodium..ppb., type = "l")
# plot(SPC14_daily$Date, SPC14_daily$Sodium..ppb., type = "l")
# 
# plot(SPC14$Date, SPC14$Magnesium..ppb., type = "l")
# plot(SPC14_daily$Date, SPC14_daily$Magnesium..ppb., type = "l")
# 
# # this is pretty questionable
# plot(SPC14$Date, SPC14$Calcium..ppb., type = "l")
# plot(SPC14_daily$Date, SPC14_daily$Calcium..ppb., type = "l")


