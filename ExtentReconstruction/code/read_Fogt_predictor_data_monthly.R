# read Fogt monthly temperature data
# read from netcdf files
library(ncdf4)
library(dplyr)

# read all data
files_mon_full_2020 = list.files("../data/fogt_predictor_data/predictor_data_through_2020", 
                            full.names = TRUE, pattern = ".nc")
files_mon_full_2021 = list.files("../data/fogt_predictor_data/updated_may2021_mon_predictor_data", 
                            full.names = TRUE, pattern = ".nc")

files_mon_2020 = list.files("../data/fogt_predictor_data/predictor_data_through_2020", 
                       full.names = FALSE, pattern = ".nc")
files_mon_2021 = list.files("../data/fogt_predictor_data/updated_may2021_mon_predictor_data", 
                           full.names = FALSE, pattern = ".nc")


#what is the difference between 2020 and 2021 data
files_mon_2021[!(files_mon_2021 %in% files_mon_2020)]
files_mon_2020[!(files_mon_2020 %in% files_mon_2021)]
# the naming convention for SSTs Nino seems to have changed


# only take data that is not in update
files_mon_full = files_mon_full_2020[!(files_mon_2020 %in% files_mon_2021)]
# manually take out SSTsNino
files_mon_full = files_mon_full[-c(97:99)]
files_mon_full = c(files_mon_full, files_mon_full_2021)

files_mon = files_mon_2020[!(files_mon_2020 %in% files_mon_2021)]
# manually take out SSTsNino
files_mon = files_mon[-c(97:99)]
files_mon = c(files_mon, files_mon_2021)


files_mon_stations = strsplit(files_mon, c(".nc")) %>% unlist()

data_list_mon = list()
# i = 6
# 162...SST NINO are problematic
for (i in c(1:5, 7:161)) {
  this.file = files_mon_full[i]
  
  # this.nc = nc_open(this.file)
  # this.var = ncvar_get(this.nc)
  
  this.raster = raster::raster(this.file)
  this.df = raster::as.data.frame(this.raster, xy = TRUE)
  names(this.df) = c("Month", "Year", files_mon_stations[i])
  data_list_mon[[i]] = this.df
}
i = 6
this.file = files_mon_full[i]
this.nc = nc_open(this.file)
this.var = ncvar_get(this.nc)

this.df = data.frame(Month = rep(1:12, 71),
                     Year = rep(1951:2021, each = 12),
                     var = c(this.var))
names(this.df)[3] = files_mon_stations[6]
data_list_mon[[6]] = this.df

# the SST Nino files are different
for (i in 162:165) {
  # print(i)
  this.file = files_mon_full[i]
  this.nc = nc_open(this.file)
  this.var = ncvar_get(this.nc)
  this.Date = as.Date("1800-1-1") + this.nc$dim$sst$vals
  
  this.df = data.frame(Month = substr(this.Date, 6, 7) %>% as.numeric(),
                       Year = substr(this.Date, 1, 4) %>% as.numeric(),
                       var = c(this.var))
  names(this.df)[3] = files_mon_stations[i]
  data_list_mon[[i]] = this.df
}


# data_list_mon = lapply(files_mon_full, 
#                            read.table, header = FALSE, fill = FALSE, 
#                            na.strings = c("-999.9", "-999.99", "-99.99", "99.99", "-9999.0", "9999.9", "7-999.92787"))
# 
# # create wide transformation of data set
# # create list of data frames
# # add station name to column names of measurement vars
# for (i in 1:length(data_list_mon)) {
#   names(data_list_mon[[i]]) = c("Year", 1:12)
#   # get measure vars
#   measure_vars = 2:13
#   this_station = strsplit(files_mon[i], ".txt")[[1]]
#   this.mdat = reshape2::melt(data_list_mon[[i]], id.vars = "Year")
#   names(this.mdat) = c("Year", "Month", this_station)
#   this.mdat = this.mdat[order(this.mdat$Year, this.mdat$Month),]
#   data_list_mon[[i]] = this.mdat
# }

# create data frame
data_mon = data_list_mon[[1]]
for (i in 2:length(data_list_mon)) {
  data_mon = full_join(data_mon, data_list_mon[[i]], by = c("Year", "Month"))
}
names(data_mon)


# compute date as mid month
data_mon$Date = apply(data_mon[, c("Year", "Month")], 1, paste0, collapse = "-") %>%
  paste0("-15") %>%
  as.Date(format = "%Y-%m-%d")
data_mon = data_mon %>%
  mutate(Doy = lubridate::yday(Date))

# reorder by date
data_mon = data_mon[order(data_mon$Date),]

# ##############################################################################
# # take out 857660_tmp, it is unstable over time
# plot(data_mon$'857660_tmp')
# data_mon = data_mon[, !names(data_mon) %in% c("857660_tmp")]




# get all variables that are measurements
measure_vars = ! (names(data_mon) %in% c("Year", "Month", "Day", "Date", "Doy"))



# ##############################################################################
# plot the time series

# subset to relevant time period
data_mon = data_mon[between(data_mon$Date,
                            as.Date("1899-01-01"), as.Date("2021-01-01")),]
data_mon$time = data_mon$Year + data_mon$Month / 12

# station 619960 is incorrect
# read actual St Helena data from 
# Feistel & Hagen (2003) Climatic changes in the subtropical Southeast Atlantic: the St. Helena Island Climate Index (1893-1999)
dat = read.delim("../data/fogt_predictor_data/StHelena/imputedStH.dat", sep = ",")
data_mon$'619960_tmp'[between(data_mon$Year, 1899, 2020)] = dat$impStH[between(dat$Year, 1899, 2020)]

dat = read.csv("~/RECONSTRUCTION/data/fogt_predictor_data/StHelena/HIX_data_no_header.dat", sep="")
head(dat)
data_mon$'619960_SLP'[between(data_mon$Year, 1899, 2020)] = dat$SLP[between(dat$Year, 1899, 2020)]

# # do an adjustment for St Helena
# raw_619010_tmp = data_mon$'619010_tmp'
# # loc_ind = data_mon$Date > as.Date("1985-01-01")
# loc_ind = (data_mon$Date > as.Date("1976-09-01"))
# loc_ind2 = loc_ind + (data_mon$Date > as.Date("1985-01-01"))
# tapply(raw_619010_tmp, loc_ind2, mean, na.rm = TRUE)
# tapply(raw_619010_tmp, loc_ind2, sd, na.rm = TRUE)
# 
# plot(data_mon$Date, raw_619010_tmp, col = loc_ind + 1)
# # create linear model for adjustment
# mod_data = cbind(data_mon, 
#                  loc_ind = factor(loc_ind + 1),
#                  loc_ind2 = factor(loc_ind2 + 1),
#                  raw_619010_tmp = raw_619010_tmp)
# limo_619010_tmp = mgcv::gam(raw_619010_tmp ~ s(time, k = 4) + 
#                               loc_ind + s(Month, by = loc_ind, bs = 'cc'), 
#                       data = mod_data)
# summary(limo_619010_tmp)
# plot(limo_619010_tmp)
# plot(data_mon$Date, mod_data$raw_619010_tmp - predict(limo_619010_tmp, newdata = mod_data))
# tapply(predict(limo_619010_tmp, newdata = mod_data) - mod_data$raw_619010_tmp,
#        loc_ind + 1, 
#        sd, na.rm = TRUE)
# # is the shift upwards observed in 1970-1980ish credible?
# plot(data_mon$Date[between(data_mon$Date,
#                            as.Date("1965-01-01"), as.Date("1982-01-01"))],
#      data_mon$'619010_tmp'[between(data_mon$Date,
#                                    as.Date("1965-01-01"), as.Date("1982-01-01"))],
#      type = "l", xlab = "Date", ylab = "St Helena Temp")

# # location change looks reasonable
# plot(data_mon$time,
#      data_mon$'619010_tmp',
#      type = "l", xlab = "Date", ylab = "St Helena Temp", xlim = c(1960, 1980))
# lines(dat$DecYear, dat$impStH, type = "l", col = 2)
# # lets look at interpolation
# plot(data_mon$time,
#      data_mon$'619010_tmp',
#      type = "l", xlab = "Date", ylab = "St Helena Temp", xlim = c(1980, 1990))
# lines(dat$DecYear, dat$impStH, type = "l", col = 2)

# # check nearby island Ascencion
# a6 <- read.csv("../data/fogt_predictor_data/Ascencion/ds570.0_monthly-6.csv")
# year.6 <- as.numeric(substr(a6[,2],1,4))
# month.6 <- as.numeric(substr(a6[,2],6,7))
# tdate.m.6 <- year.6 + month.6/12 - 0.5/12
# plot(x=tdate.m.6[year.6 < 1977],a6[year.6 < 1977,"T.degC."], xlim = c(1965, 1982), type = "l")
# # do a monthly detrend on this data
# a6$Month = month.6
# limo_a6 = mgcv::gam(T.degC. ~ s(Month, bs = 'cc'), 
#                             data = a6[year.6 < 1977,])
# plot(limo_a6)
# summary(limo_a6)
# plot(x=tdate.m.6[year.6 < 1977],
#      a6[year.6 < 1977,"T.degC."] - predict(limo_a6, newdata = a6[year.6 < 1977,]),
#      xlim = c(1965, 1982), type = "l",
#      xlab = 'Year', ylab = "Ascencion Temp - Monthly Avg Temp")



# 837430_tmp has an implausible strong decline in last couple years?
# 919380_tmp seems to have too much of a trend
for (i in 1:sum(measure_vars)) { # go from 1 to 
  this.var = names(data_mon)[measure_vars][i]
  first_non_miss = which(!is.na(data_mon[,this.var]))[1]
  # png(paste0("../plots/covariates/raw_covariates/", this.var, ".jpg"), width = 720, height = 480)
  plot(data_mon$Date[first_non_miss:nrow(data_mon)],
       data_mon[first_non_miss:nrow(data_mon),this.var],
       type = "l", main = this.var)
  # dev.off()
}

# take out implausible data
data_mon$'837430_tmp'[data_mon$Date > as.Date("2017-05-01")] = NA
data_mon = data_mon[, !(names(data_mon) %in% c("919380_tmp"))]
# recompute which variables are measurements
measure_vars = ! (names(data_mon) %in% c("Year", "Month", "Day", "Date", "Doy"))


# ##############################################################################
# compute missingness
# per day
miss_measurements_per_day = apply(data_mon[,measure_vars], 1, function(x) {
  mean(is.na(x)) * 100
})
# png("../plots/covariates/miss_per_day.png", width = 720, height = 480)
hist(miss_measurements_per_day, main = "% missing measurements per day", xlab = "% missing")
dev.off()
# per station
miss_measurements_per_station = apply(data_mon[,measure_vars], 2, function(x) {
  mean(is.na(x)) * 100
})

# png("../plots/covariates/miss_per_station.png", width = 720, height = 480)
hist(miss_measurements_per_station, main = "% missing measurements per station", xlab = "% missing")
dev.off()

# missing measurement per day and station
miss_meas = apply(data_mon[,measure_vars], 2, is.na)
# png("../plots/covariates/missing_mon.png", width = 720, height = 480)
image(x = data_mon$Date, y = 1:ncol(miss_meas), z = miss_meas, 
      xlab = "Day", ylab = "Station", 
      xlim = c(as.Date("1900-01-01"), as.Date("2021-06-01")),
      col = c("black", "white"))
dev.off()

# ##############################################################################
# get the raw correlations
raw_cor = cor(data_mon[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/raw_correlations_mon.png", width = 2000, height = 2000)
corrplot::corrplot(raw_cor, method="color")
dev.off()

# detrend data
# get all variables that are measurements
measure_vars = ! (names(data_mon) %in% c("Year", "Month", "Day", "Date", "Doy"))
data_mon_detrended = data_mon
Vdata_mon = data_mon
names(Vdata_mon) = paste0("V", names(data_mon))
# use this instead
# gam(ts.nsidc ~ tdate.m + s(ax.nophase, bs="cc",k=6,fx=FALSE),knots=list(ax=c(0,12)))
for (i in 1:sum(measure_vars)) {
  this.var = names(data_mon)[measure_vars][i]
  # this.gam = mgcv::gam(formula(paste0("V", this.var, " ~ s(VMonth, bs = 'cc')")),
  #                      data = Vdata_mon)
  this.gam = mgcv::gam(formula(paste0("V", this.var, " ~ s(VMonth, bs = 'cc', k = 6)")),
                       knots = list(VMonth = seq(0.5, 12.5, length.out = 6)),
                       data = Vdata_mon)
  data_mon_detrended[,this.var] = data_mon[,this.var] - unlist(predict(this.gam, newdata = Vdata_mon))
  data_mon_detrended[,this.var] = c(scale(data_mon_detrended[,this.var], center = TRUE, scale = TRUE))
}

# get the detrended correlations
detrended_cor = cor(data_mon_detrended[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/detrended_correlations_mon.png", width = 2000, height = 2000)
corrplot::corrplot(detrended_cor, method="color")
dev.off()

