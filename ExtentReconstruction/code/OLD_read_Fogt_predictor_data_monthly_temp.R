# read Fogt monthly temperature data
library(dplyr)

# delete first two lines of Byrd station data, 891250_tmp
# read all data
files_monTemp_full = list.files("../data/fogt_predictor_data/data_txt_files/mon_txt_tmp", 
                                full.names = TRUE, pattern = ".txt")
files_monTemp = list.files("../data/fogt_predictor_data/data_txt_files/mon_txt_tmp", 
                           full.names = FALSE, pattern = ".txt")
files_monTemp_full = files_monTemp_full[-c(58, 89)]
files_monTemp = files_monTemp[-c(58, 89)]

# delete station_temp_data.txt
# I am not really sure what this is anyways
# delete 916430_tmp, this seems unreliable alltogether

files_monTemp_stations = strsplit(files_monTemp, ".txt")
data_list_monTemp = lapply(files_monTemp_full, 
                           read.table, header = FALSE, fill = FALSE, 
                           na.strings = c("-999.9", "-999.99", "-99.99", "99.99", "-9999.0", "9999.9", "7-999.92787"))


# create wide transformation of data set
# create list of data frames
# add station name to column names of measurement vars
for (i in 1:length(data_list_monTemp)) {
  names(data_list_monTemp[[i]]) = c("Year", 1:12)
  # get measure vars
  measure_vars = 2:13
  this_station = strsplit(files_monTemp[i], ".txt")[[1]]
  this.mdat = reshape2::melt(data_list_monTemp[[i]], id.vars = "Year")
  names(this.mdat) = c("Year", "Month", this_station)
  this.mdat = this.mdat[order(this.mdat$Year, this.mdat$Month),]
  data_list_monTemp[[i]] = this.mdat
}

# create data frame
data_monTemp = data_list_monTemp[[1]]
for (i in 2:length(files_monTemp)) {
  data_monTemp = full_join(data_monTemp, data_list_monTemp[[i]], by = c("Year", "Month"))
}
names(data_monTemp)

# compute day of year
# delete February 29, 1900, which is a day that does not exist
data_monTemp$Date = apply(data_monTemp[, c("Year", "Month")], 1, paste0, collapse = "-") %>%
  paste0("-15") %>%
  as.Date(format = "%Y-%m-%d")
data_monTemp = data_monTemp %>%
  mutate(Doy = lubridate::yday(Date))
# get all variables that are measurements
measure_vars = ! (names(data_monTemp) %in% c("Year", "Month", "Day", "Date", "Doy"))
# reorder by date
data_monTemp = data_monTemp[order(data_monTemp$Date),]

# ##############################################################################
# plot the time series

data_monTemp$'946370_tmp'[data_monTemp$'946370_tmp' < 1] = NA
data_monTemp$'944300_tmp'[data_monTemp$'944300_tmp' < 1] = NA
data_monTemp$'943460_tmp'[data_monTemp$'943460_tmp' < 1] = NA
data_monTemp$'943120_tmp'[data_monTemp$'943120_tmp' < 1] = NA
data_monTemp$'943000_tmp'[data_monTemp$'943000_tmp' < 1] = NA
data_monTemp$'943000_tmp'[data_monTemp$'943000_tmp' < 1] = NA
data_monTemp$'939870_tmp'[data_monTemp$'939870_tmp' < 1] = NA
data_monTemp$'938940_tmp'[data_monTemp$'938940_tmp' < 1] = NA
data_monTemp$'931190_tmp'[data_monTemp$'931190_tmp' < 5] = NA
data_monTemp$'916800_tmp'[data_monTemp$'916800_tmp' < 1] = NA
# # take out 916430_tmp data seems unreliable
# data_monTemp$'916430_tmp'[data_monTemp$'916430_tmp' < 25 |
#                             data_monTemp$'916430_tmp' > 35] = NA
data_monTemp$'873440_tmp'[data_monTemp$'873440_tmp' < 1] = NA
# 857660_tmp not plausible, maybe cut in half for before and after 1940?
data_monTemp$'857660_tmp'[data_monTemp$'Year' <= 1941] = NA
data_monTemp$'854420_tmp'[data_monTemp$'854420_tmp' < 1] = NA
data_monTemp$'837430_tmp'[data_monTemp$'837430_tmp' < 15] = NA
data_monTemp$'688580_tmp'[data_monTemp$'688580_tmp' < 1] = NA
data_monTemp$'685880_tmp'[data_monTemp$'685880_tmp' < 1] = NA
data_monTemp$'688160_tmp'[data_monTemp$'688160_tmp' < 1] = NA
data_monTemp$'619980_tmp'[data_monTemp$'619980_tmp' > 15] = NA
# 619010_tmp seems implausible, before and after gap in 1980s the values changedß
data_monTemp$'619010_tmp'[data_monTemp$'619010_tmp' > 22] = NA

for (i in 1:sum(measure_vars)) {
  this.var = names(data_monTemp)[measure_vars][i]
  first_non_miss = which(!is.na(data_monTemp[,this.var]))[1]
  plot(data_monTemp$Date[first_non_miss:nrow(data_monTemp)],
       data_monTemp[first_non_miss:nrow(data_monTemp),this.var],
       type = "l", main = this.var)
}

# ##############################################################################
# compute missingness
# per day
miss_measurements_per_day = apply(data_monTemp[,measure_vars], 1, function(x) {
  mean(is.na(x)) * 100
})
# png("../plots/covariates/miss_per_day.png", width = 720, height = 480)
hist(miss_measurements_per_day, main = "% missing measurements per day", xlab = "% missing")
dev.off()
# per station
miss_measurements_per_station = apply(data_monTemp[,measure_vars], 2, function(x) {
  mean(is.na(x)) * 100
})

# png("../plots/covariates/miss_per_station.png", width = 720, height = 480)
hist(miss_measurements_per_station, main = "% missing measurements per day", xlab = "% missing")
dev.off()

# missing measurement per day and station
miss_meas = apply(data_monTemp[,measure_vars], 2, is.na)
# png("../plots/covariates/missing_monTemp.png", width = 720, height = 480)
image(x = data_monTemp$Date, y = 1:ncol(miss_meas), z = miss_meas, xlab = "Day", ylab = "Station", 
      col = c("black", "white"))
dev.off()

# ##############################################################################
# get the raw correlations
raw_cor = cor(data_monTemp[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/raw_correlations_monTemp.png", width = 2000, height = 2000)
corrplot::corrplot(raw_cor, method="color")
dev.off()

# detrend data
data_monTemp_detrended = data_monTemp
Vdata_monTemp = data_monTemp
names(Vdata_monTemp) = paste0("V", names(data_monTemp))
for (i in 1:sum(measure_vars)) {
  this.var = names(data_monTemp)[measure_vars][i]
  this.gam = mgcv::gam(formula(paste0("V", this.var, " ~ s(VDoy, bs = 'cc')")), data = Vdata_monTemp)
  data_monTemp_detrended[,this.var] = data_monTemp[,this.var] - predict(this.gam, newdata = Vdata_monTemp)
  data_monTemp_detrended[,this.var] = scale(data_monTemp_detrended[,this.var], center = TRUE, scale = TRUE)
}

# get the detrended correlations
detrended_cor = cor(data_monTemp_detrended[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/detrended_correlations_monTemp.png", width = 2000, height = 2000)
corrplot::corrplot(detrended_cor, method="color")
dev.off()