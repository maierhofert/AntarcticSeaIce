# read Fogt daily data
# read txt files
# why is Feb 29, 1900 in the data set?
# all values are missing but it shouldn't exist nonetheless

# changed ASN00034002 column name MAX to TMAX

library(dplyr)

# read all data
files_NCDC_full = list.files("../data/fogt_predictor_data/data_txt_files/daily_txt_NCDC_tmp_slp", full.names = TRUE)
files_NCDC = list.files("../data/fogt_predictor_data/data_txt_files/daily_txt_NCDC_tmp_slp", full.names = FALSE)
files_NCDC_stations = strsplit(files_NCDC, ".txt")
data_list_NCDC = lapply(files_NCDC_full, 
                         read.table, header = FALSE, fill = TRUE, 
                        na.strings = c("-999.9", "-9999.0", "9999.9", "7-999.92787"))

# create wide transformation of data set
# create list of data frames
# add station name to column names of measurement vars
for (i in 1:length(files_NCDC)) {
  names(data_list_NCDC[[i]]) = c("Year", "Month", "Day", "Temp", "SLP")
  # get measure vars
  measure_vars = which(names(data_list_NCDC[[i]]) %in% c("Temp", "SLP"))
  this_station = strsplit(files_NCDC[i], ".txt")[[1]]
  names(data_list_NCDC[[i]])[measure_vars] = paste0(names(data_list_NCDC[[i]])[measure_vars], "_", this_station)
}
# there are some coding errors in here
data_list_NCDC[[10]]$Temp_837430 = as.numeric(data_list_NCDC[[10]]$Temp_837430) 
# which(is.na(as.numeric(data_list_NCDC[[10]]$Temp_837430)) & !is.na(data_list_NCDC[[10]]$Temp_837430))


# create data frame
data_NCDC = data_list_NCDC[[1]]
for (i in 2:length(files_NCDC)) {
  data_NCDC = full_join(data_NCDC, data_list_NCDC[[i]], by = c("Year", "Month", "Day"))
}
names(data_NCDC)

# compute day of year
# delete February 29, 1900, which is a day that does not exist
data_NCDC$Date = apply(data_NCDC[, c("Year", "Month", "Day")], 1, paste0, collapse = "-") %>% as.Date()
data_NCDC = data_NCDC %>%
  mutate(Doy = lubridate::yday(Date)) %>%
  filter(!(Day == 29 & Month == 2 & Year == 1900))
# SLP_893240 entirely missing
data_NCDC = select(data_NCDC, !SLP_893240)
# get all variables that are measurements
measure_vars = ! (names(data_NCDC) %in% c("Year", "Month", "Day", "Date", "Doy"))
# reorder by date
data_NCDC = data_NCDC[order(data_NCDC$Date),]

# ##############################################################################
# plot the time series
# SLP look often truncated
# SLP_890550, SLP_855850,  extreme values to the top look very weird, are all one of two values
data_NCDC$SLP_837430[data_NCDC$SLP_837430 >= 10000] = NA
data_NCDC$Temp_837430[data_NCDC$Temp_837430 <= 0] = NA
# SLP_619010 looks crazy, Temp_619010 looks implausible as well

for (i in 1:sum(measure_vars)) {
  this.var = names(data_NCDC)[measure_vars][i]
  first_non_miss = which(!is.na(data_NCDC[,this.var]))[1]
  plot(data_NCDC$Date[first_non_miss:nrow(data_NCDC)],
       data_NCDC[first_non_miss:nrow(data_NCDC),this.var],
       type = "l", main = this.var)
}

# ##############################################################################
# compute missingness
# per day
miss_measurements_per_day = apply(data_NCDC[,measure_vars], 1, function(x) {
  mean(is.na(x)) * 100
})
# png("../plots/covariates/miss_per_day.png", width = 720, height = 480)
hist(miss_measurements_per_day, main = "% missing measurements per day", xlab = "% missing")
dev.off()
# per station
miss_measurements_per_station = apply(data_NCDC[,measure_vars], 2, function(x) {
  mean(is.na(x)) * 100
})

# png("../plots/covariates/miss_per_station.png", width = 720, height = 480)
hist(miss_measurements_per_station, main = "% missing measurements per day", xlab = "% missing")
dev.off()

# missing measurement per day and station
miss_meas = apply(data_NCDC[,measure_vars], 2, is.na)
# png("../plots/covariates/missing_NCDC.png", width = 720, height = 480)
image(x = data_NCDC$Date, y = 1:ncol(miss_meas), z = miss_meas, xlab = "Day", ylab = "Station", 
      col = c("black", "white"))
dev.off()

# ##############################################################################
# get the raw correlations
raw_cor = cor(data_NCDC[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/raw_correlations_SLP.png", width = 2000, height = 2000)
corrplot::corrplot(raw_cor, method="color")
dev.off()

# detrend data
data_NCDC_detrended = data_NCDC
for (i in 1:sum(measure_vars)) {
  this.var = names(data_NCDC)[measure_vars][i]
  this.gam = mgcv::gam(formula(paste0(this.var, " ~ s(Doy, bs = 'cc')")), data = data_NCDC)
  data_NCDC_detrended[,this.var] = data_NCDC[,this.var] - predict(this.gam, newdata = data_NCDC)
  data_NCDC_detrended[,this.var] = scale(data_NCDC_detrended[,this.var], center = TRUE, scale = TRUE)
}

# get the detrended correlations
detrended_cor = cor(data_NCDC_detrended[,measure_vars], use = "pairwise.complete.obs")
# png("../plots/covariates/detrended_correlations_SLP.png", width = 2000, height = 2000)
corrplot::corrplot(detrended_cor, method="color")
dev.off()