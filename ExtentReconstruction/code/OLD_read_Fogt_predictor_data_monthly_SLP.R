# read Fogt monthly data
library(dplyr)

# read all data
files_monSLP_full = list.files("../data/fogt_predictor_data/data_txt_files/mon_txt_SLP", 
                               full.names = TRUE, pattern = ".txt")
files_monSLP = list.files("../data/fogt_predictor_data/data_txt_files/mon_txt_SLP", 
                          full.names = FALSE, pattern = ".txt")
files_monSLP_stations = strsplit(files_monSLP, ".txt")
data_list_monSLP = lapply(files_monSLP_full, 
                         read.table, header = FALSE, fill = FALSE, 
                        na.strings = c("-999.9", "-9999.0", "9999.9", "7-999.92787"))
# create wide transformation of data set
# create list of data frames
# add station name to column names of measurement vars
for (i in 1:length(files_monSLP)) {
  names(data_list_monSLP[[i]]) = c("Year", 1:12)
  # get measure vars
  measure_vars = 2:13
  this_station = strsplit(files_monSLP[i], ".txt")[[1]]
  this.mdat = reshape2::melt(data_list_monSLP[[i]], id.vars = "Year")
  names(this.mdat) = c("Year", "Month", this_station)
  this.mdat = this.mdat[order(this.mdat$Year, this.mdat$Month),]
  data_list_monSLP[[i]] = this.mdat
}

# create data frame
data_monSLP = data_list_monSLP[[1]]
for (i in 2:length(files_monSLP)) {
  data_monSLP = full_join(data_monSLP, data_list_monSLP[[i]], by = c("Year", "Month"))
}
names(data_monSLP)

# compute day of year
# delete February 29, 1900, which is a day that does not exist
data_monSLP$Date = apply(data_monSLP[, c("Year", "Month")], 1, paste0, collapse = "-") %>%
  paste0("-15") %>%
  as.Date(format = "%Y-%m-%d")
data_monSLP = data_monSLP %>%
  mutate(Doy = lubridate::yday(Date))
# get all variables that are measurements
measure_vars = ! (names(data_monSLP) %in% c("Year", "Month", "Day", "Date", "Doy"))
# reorder by date
data_monSLP = data_monSLP[order(data_monSLP$Date),]

# ##############################################################################
# plot the time series
# IS 954880 plausible?
# Scott base McMurdo SLP is implausible between 2013 and 2017
data_monSLP$Scott_Base_McMurdo_SLP[between(data_monSLP$Date, 
                                           as.Date("2013-08-01"), as.Date("2017-01-01"))] = NA
data_monSLP$'943260_SLP'[data_monSLP$'943260_SLP' > 2000] = NA
data_monSLP$'916800_SLP'[data_monSLP$'916800_SLP' < 100 |
                           data_monSLP$'916800_SLP' > 1019] = NA
data_monSLP$'891250_SLP'[data_monSLP$'891250_SLP' < 770] = NA
data_monSLP$'890550_SLP'[data_monSLP$'890550_SLP' < 970] = NA
data_monSLP$'889630_SLP'[data_monSLP$'889630_SLP' > 9796] = NA
data_monSLP$'879380_SLP'[data_monSLP$'879380_SLP' < 980] = NA
data_monSLP$'878490_SLP'[data_monSLP$'878490_SLP' < 110] = NA
data_monSLP$'875440_SLP'[data_monSLP$'875440_SLP' > 1030] = NA
data_monSLP$'872220_SLP'[data_monSLP$'872220_SLP' < 200] = NA
data_monSLP$'870470_SLP'[data_monSLP$'870470_SLP' > 900] = NA
data_monSLP$'688580_tmp'[data_monSLP$'688580_tmp' < 10] = NA
data_monSLP$'688580_SLP'[data_monSLP$'688580_SLP' < 200 |
                           data_monSLP$'688580_SLP' < 1011.1] = NA
data_monSLP$'688160_SLP'[data_monSLP$'688160_SLP' < 0] = NA
data_monSLP$'685880_SLP'[data_monSLP$'685880_SLP' < 200] = NA

for (i in 1:sum(measure_vars)) {
  this.var = names(data_monSLP)[measure_vars][i]
  first_non_miss = which(!is.na(data_monSLP[,this.var]))[1]
  plot(data_monSLP$Date[first_non_miss:nrow(data_monSLP)],
       data_monSLP[first_non_miss:nrow(data_monSLP),this.var],
       type = "l", main = this.var)
}

# ##############################################################################
# compute missingness
# per day
miss_measurements_per_day = apply(data_monSLP[,measure_vars], 1, function(x) {
  mean(is.na(x)) * 100
})
# png("../plots/covariates/miss_per_day.png", width = 720, height = 480)
hist(miss_measurements_per_day, main = "% missing measurements per day", xlab = "% missing")
dev.off()
# per station
miss_measurements_per_station = apply(data_monSLP[,measure_vars], 2, function(x) {
  mean(is.na(x)) * 100
})

# png("../plots/covariates/miss_per_station.png", width = 720, height = 480)
hist(miss_measurements_per_station, main = "% missing measurements per day", xlab = "% missing")
dev.off()

# missing measurement per day and station
miss_meas = apply(data_monSLP[,measure_vars], 2, is.na)
# png("../plots/covariates/missing_monSLP.png", width = 720, height = 480)
image(x = data_monSLP$Date, y = 1:ncol(miss_meas), z = miss_meas, xlab = "Day", ylab = "Station", 
      col = c("black", "white"))
dev.off()

# ##############################################################################
# get the raw correlations
raw_cor = cor(data_monSLP[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/raw_correlations_monSLP.png", width = 2000, height = 2000)
corrplot::corrplot(raw_cor, method="color")
dev.off()

# detrend data
data_monSLP_detrended = data_monSLP
Vdata_monSLP = data_monSLP
names(Vdata_monSLP) = paste0("V", names(data_monSLP))
for (i in 1:sum(measure_vars)) {
  this.var = names(data_monSLP)[measure_vars][i]
  this.gam = mgcv::gam(formula(paste0("V", this.var, " ~ s(VDoy, bs = 'cc')")), data = Vdata_monSLP)
  data_monSLP_detrended[,this.var] = data_monSLP[,this.var] - predict(this.gam, newdata = Vdata_monSLP)
  data_monSLP_detrended[,this.var] = scale(data_monSLP_detrended[,this.var], center = TRUE, scale = TRUE)
}

# get the detrended correlations
detrended_cor = cor(data_monSLP_detrended[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/detrended_correlations_monSLP.png", width = 2000, height = 2000)
corrplot::corrplot(detrended_cor, method="color")
dev.off()