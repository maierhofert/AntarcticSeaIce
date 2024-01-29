# read Fogt daily data
# read txt files
# why is Feb 29, 1900 in the data set?
# all values are missing but it shouldn't exist nonetheless

# changed ASN00034002 column name MAX to TMAX

library(dplyr)

# # read daily_txt_ghcnd_tmp
# datASN2011 = read.table("../data/fogt_predictor_data/data_txt_files/daily_txt_ghcnd_tmp/ASN00002011.txt", 
#                         header = TRUE, na.strings = c("-999.9", "-9999.0"))
# datASN2011$Date = apply(datASN2011[, c("Year", "Month", "Day")], 1, paste0, collapse = "-") %>% as.Date()
# # compute day of year
# # delete February 29, 1900, which is a day that does not exist
# datASN2011 = datASN2011 %>%
#   mutate(Doy = lubridate::yday(Date)) %>%
#   filter(!(Day == 29 & Month == 2 & Year == 1900))
# datASN2011
# 
# head(datASN2011)
# tail(datASN2011)
# # plot data
# plot(datASN2011$Date, datASN2011$TMIN, type = "l")
# # we have strong seasonality, some missing data
# plot(datASN2011$Date[1:1000], datASN2011$TMIN[1:1000], type = "l", ylim = c(30, 110))
# lines(datASN2011$Date[1:1000], datASN2011$TMAX[1:1000], type = "l", col = 2)
# 
# 
# # some days have missing values in only tmin or tmax, this could be used for filling NAs
# table(is.na(datASN2011$TMAX), is.na(datASN2011$TMIN))
# # strong correlation between daily min and max temperature
# cor(datASN2011$TMAX, datASN2011$TMIN, use = "pairwise.complete.obs") %>% round(2)
# 
# # some other data
# datASN3002 = read.table("data/fogt_predictor_data/data_txt_files/daily_txt_ghcnd_tmp/ASN00003002.txt", 
#                         header = TRUE, na.strings = "-999.9")
# head(datASN3002)
# tail(datASN3002)
# summary(datASN3002)
# 
# 
# # next data set
# datASN4020 = read.table("data/fogt_predictor_data/data_txt_files/daily_txt_ghcnd_tmp/ASN00003002.txt", 
#                         header = TRUE, na.strings = "-999.9")
# head(datASN4020)
# tail(datASN4020)
# summary(datASN4020)


# read all data
files_ghcnd_full = list.files("../data/fogt_predictor_data/data_txt_files/daily_txt_ghcnd_tmp", full.names = TRUE)
files_ghcnd = list.files("../data/fogt_predictor_data/data_txt_files/daily_txt_ghcnd_tmp", full.names = FALSE)
files_ghcnd
files_ghcnd_stations = strsplit(files_ghcnd, ".txt")
data_list_ghcnd = lapply(files_ghcnd_full, 
                         read.table, header = TRUE, na.strings = c("-999.9", "-9999.0"))

# create wide transformation of data set
# create list of data frames
# add station name to column names of measurement vars
for (i in 1:length(files_ghcnd)) {
  names(data_list_ghcnd[[i]])[names(data_list_ghcnd[[i]]) %in% c("Mon", "mon")] = "Month"
  names(data_list_ghcnd[[i]])[names(data_list_ghcnd[[i]]) == "year"] = "Year"
  # get measure vars
  measure_vars = which(names(data_list_ghcnd[[i]]) %in% c("TAVG", "TMAX", "TMIN"))
  this_station = strsplit(files_ghcnd[i], ".txt")[[1]]
  names(data_list_ghcnd[[i]])[measure_vars] = paste0(names(data_list_ghcnd[[i]])[measure_vars], "_", this_station)
}

# create data frame
data_ghcnd = data_list_ghcnd[[1]]
for (i in 2:length(files_ghcnd)) {
  data_ghcnd = full_join(data_ghcnd, data_list_ghcnd[[i]], by = c("Year", "Month", "Day"))
}
names(data_ghcnd)

# compute day of year
# delete February 29, 1900, which is a day that does not exist
data_ghcnd$Date = apply(data_ghcnd[, c("Year", "Month", "Day")], 1, paste0, collapse = "-") %>% as.Date()
data_ghcnd = data_ghcnd %>%
  mutate(Doy = lubridate::yday(Date)) %>%
  filter(!(Day == 29 & Month == 2 & Year == 1900))
measure_vars = ! (names(data_ghcnd) %in% c("Year", "Month", "Day", "Date", "Doy"))

# ##############################################################################
# plot the time series
# there are problems with outliers
data_ghcnd$TMAX_ASN00028004[data_ghcnd$TMAX_ASN00028004 < 40] = NA 
# are some missing values set to exactly zero, e.g. in ASN00074128?
data_ghcnd$TMIN_ASN00074128[data_ghcnd$TMIN_ASN00074128 < 1 | is.na(data_ghcnd$TMIN_ASN00074128)] = NA 
# TMIN_ASN00040264 has some unplausible looking years around the gap in the 1930s
# set 0.0 values in ASN00074128 as missing
data_ghcnd$TMIN_ASN00074128[data_ghcnd$TMIN_ASN00074128 < 1 | is.na(data_ghcnd$TMIN_ASN00074128)] = NA 
data_ghcnd$TMAX_ASN00074128[data_ghcnd$TMAX_ASN00074128 < 1 | is.na(data_ghcnd$TMAX_ASN00074128)] = NA 

data_ghcnd$TMIN_AYW00090001[data_ghcnd$TMIN_AYW00090001 > 10] = NA 
data_ghcnd$TMAX_AYW00090001[data_ghcnd$TMAX_AYW00090001 > 10] = NA
data_ghcnd$TMIN_AYM00088963[data_ghcnd$TMIN_AYM00088963 > 50] = NA 
# TMIN_ASN00040264 strange hole in measurements with implausible values around
data_ghcnd$TMIN_ASN00040264[between(data_ghcnd$Date, as.Date("1924-01-01"), as.Date("1936-01-01"))] = NA 

for (i in 1:sum(measure_vars)) {
  this.var = names(data_ghcnd)[measure_vars][i]
  first_non_miss = which(!is.na(data_ghcnd[,this.var]))[1]
  plot(data_ghcnd$Date[first_non_miss:nrow(data_ghcnd)],
       data_ghcnd[first_non_miss:nrow(data_ghcnd),this.var], 
       type = "l", main = this.var)
}


# ##############################################################################
# compute missingness
# per day
miss_measurements_per_day = apply(data_ghcnd[,measure_vars], 1, function(x) {
  mean(is.na(x)) * 100
})
# png("../plots/covariates/miss_per_day.png", width = 720, height = 480)
hist(miss_measurements_per_day, main = "% missing measurements per day", xlab = "% missing")
dev.off()
# per station
miss_measurements_per_station = apply(data_ghcnd[,measure_vars], 2, function(x) {
  mean(is.na(x)) * 100
})

# png("../plots/covariates/miss_per_station.png", width = 720, height = 480)
hist(miss_measurements_per_station, main = "% missing measurements per day", xlab = "% missing")
dev.off()

# missing measurement per day and station
miss_meas = apply(data_ghcnd[,measure_vars], 2, is.na)
# png("../plots/covariates/missing_ghcnd.png", width = 720, height = 480)
image(x = data_ghcnd$Date, y = 1:ncol(miss_meas), z = miss_meas, xlab = "Day", ylab = "Station", 
      col = c("black", "white"))
dev.off()

# ##############################################################################
# get the raw correlations
raw_cor = cor(data_ghcnd[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/raw_correlations_ghcnd.png", width = 2000, height = 2000)
corrplot::corrplot(raw_cor, method="color")
dev.off()

# detrend data
data_ghcnd_detrended = data_ghcnd
for (i in 1:sum(measure_vars)) {
  this.var = names(data_ghcnd)[measure_vars][i]
  this.gam = mgcv::gam(formula(paste0(this.var, " ~ s(Doy, bs = 'cc')")), data = data_ghcnd)
  data_ghcnd_detrended[,this.var] = data_ghcnd[,this.var] - predict(this.gam, newdata = data_ghcnd)
  data_ghcnd_detrended[,this.var] = scale(data_ghcnd_detrended[,this.var], center = TRUE, scale = TRUE)
}
plot(data_ghcnd_detrended$Date, data_ghcnd_detrended$TMAX_ASN00002011, type = "l")

# get the detrended correlations
detrended_cor = cor(data_ghcnd_detrended[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/detrended_correlations_ghcnd.png", width = 2000, height = 2000)
corrplot::corrplot(detrended_cor, method="color")
dev.off()




