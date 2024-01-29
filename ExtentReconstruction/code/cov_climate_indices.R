library(ncdf4)
# read large scale climate indices
files_mon_full = list.files("../data/fogt_predictor_data/monthly_predictor_data", 
                            full.names = TRUE, pattern = ".nc")
files_mon = list.files("../data/fogt_predictor_data/monthly_predictor_data/", 
                       full.names = FALSE, pattern = ".nc")
files_mon_stations = strsplit(files_mon, c(".nc")) %>% unlist()


data_list_mon = list()
for (i in c(1)) {
  this.file = files_mon_full[i]
  
  this.nc = nc_open(this.file)
  this.var = ncvar_get(this.nc)
  this.time = ncvar_get(this.nc, "prmsl")
  this.time = as.POSIXct(this.time*60^2, origin='1800-01-01 00:00')
  
  # plot(this.time, this.var)
  
  this.df = data.frame(Month = substr(this.time, 6, 7) %>% as.numeric(),
                       Year = substr(this.time, 1, 4) %>% as.numeric())
  this.df[,files_mon[i]] = this.var
  data_list_mon[[i]] = this.df
}
for (i in c(4, 5)) {
  this.file = files_mon_full[i]
  
  this.nc = nc_open(this.file)
  this.var = ncvar_get(this.nc)
  this.time = ncvar_get(this.nc, "msl")
  this.time = as.POSIXct(this.time*60^2, origin='1900-01-01 00:00')
  
  plot(this.time, this.var)
  this.df = data.frame(
    Month = substr(this.time, 6, 7) %>% as.numeric(),
    Year = substr(this.time, 1, 4) %>% as.numeric())
  this.df[,files_mon[i]] = this.var
  data_list_mon[[i]] = this.df
}

for (i in c(7, 11)) {
  this.file = files_mon_full[i]
  this.nc = nc_open(this.file)
  this.var = ncvar_get(this.nc)
  this.time = ncvar_get(this.nc, "yrmon")
  
  this.df = data.frame(
    Month = substr(this.time, 5, 6) %>% as.numeric(),
    Year = substr(this.time, 1, 4) %>% as.numeric())
  this.df[,files_mon[i]] = this.var
  data_list_mon[[i]] = this.df
}

for (i in c(15:18)) {
  this.file = files_mon_full[i]
  this.nc = nc_open(this.file)
  this.var = ncvar_get(this.nc)
  this.time = ncvar_get(this.nc, "sst")
  this.time = as.POSIXct(this.time*60^2*24, origin='1800-01-01 00:00')
  
  this.df = data.frame(
    Month = substr(this.time, 6, 7) %>% as.numeric(),
    Year = substr(this.time, 1, 4) %>% as.numeric())
  this.df[,files_mon[i]] = this.var
  data_list_mon[[i]] = this.df
}

for (i in c(2, 3, 6, 8, 9, 12)) {
  this.file = files_mon_full[i]
  this.raster = raster::raster(this.file)
  this.df = raster::as.data.frame(this.raster, xy = T)
  names(this.df) = c("Month", "Year", files_mon_stations[i])
  data_list_mon[[i]] = this.df
}
for (i in c(10)) {
  this.file = files_mon_full[i]
  this.raster = raster::raster(this.file)
  this.df = raster::as.data.frame(this.raster, xy = T)
  names(this.df) = c("Season", "Year", files_mon_stations[i])
  this.df$Month = this.df$Season * 3 - 2
  data_list_mon[[i]] = this.df[,-1]
}

# whats going on with 14? sst.mon.ltm.1981-2010.nc
# this is a climatology for the 12 months out of the year (seems useless)
# whats going on with 13? sst.mnmean.nc
# 13 seems to be a SST reconstruction (this could be useful?)
for (i in c(13, 14)) {
  this.file = files_mon_full[i]
  this.brick = raster::brick(this.file)
  # this.arr = 
  # image(this.arr[,,1])
  # this.df = raster::as.data.frame(this.raster, xy = T)
  # names(this.df) = c("Month", "Year", files_mon_stations[i])
  # data_list_mon[[i]] = this.df
}

str(data_list_mon)

data_list_mon = data_list_mon[-(13:14)]

res_data = data_list_mon[[1]]
for (i in 2:length(data_list_mon)) {
  res_data = full_join(res_data, data_list_mon[[i]], by = c("Month", "Year"))
}
res_data$Date = paste0(res_data$Year, "-", res_data$Month, "-15") %>% as.Date()
res_data = res_data[order(res_data$Date),]
res_data = res_data[res_data$Date > as.Date("1900-01-01"),]
str(res_data)

this.cor = cor(res_data[,-c(1, 2, 19)], use = "pairwise") %>% round(2)
heatmap(this.cor)

res_data$IPOunfilteredV5[res_data$IPOunfilteredV5 < -98] = NA
summary(res_data)

# detrend data
res_data$Doy = lubridate::yday(res_data$Date)
names(res_data)[3] = "V20cr.mslp.std.SAMIndices"
res_data_detrended = res_data
# 12 SAMRecon is seasonal, not monthly
for (i in c(3:11, 13:18)) {
  this.var = names(res_data)[i]
  this.gam = mgcv::gam(formula(paste0(this.var, " ~ s(Doy, bs = 'cc')")), data = res_data)
  res_data_detrended[,this.var] = res_data[,this.var] - unlist(predict(this.gam, newdata = res_data))
  res_data_detrended[,this.var] = c(scale(res_data_detrended[,this.var], center = TRUE, scale = TRUE))
}
# do SAMRecon
this.var = "SAMRecon"
this.gam = mgcv::gam(formula(paste0(this.var, " ~ as.factor(Month)")), data = res_data)
res_data_detrended[!is.na(res_data_detrended$SAMRecon),this.var] = res_data[!is.na(res_data_detrended$SAMRecon),this.var] - 
  unlist(predict(this.gam, newdata = res_data[!is.na(res_data_detrended$SAMRecon),]))
res_data_detrended[,this.var] = c(scale(res_data_detrended[,this.var], center = TRUE, scale = TRUE))
# this one is seasonal, not monthly, so do interpolation first
# there are *no* missing values
smooth = ksmooth(res_data_detrended$Date[!is.na(res_data_detrended$SAMRecon)], 
                 res_data_detrended[!is.na(res_data_detrended$SAMRecon),"SAMRecon"],
                 kernel = "normal", bandwidth = 30.4 * 4)$y
SplineFun <- splinefun(x = res_data_detrended$Date[!is.na(res_data_detrended$SAMRecon)],
                       y = smooth,
                       method = "natural")
res_data_detrended$SAMRecon = SplineFun(res_data_detrended$Date)
res_data_detrended$SAMRecon[1420:1440] = NA
SplineFitSAMRecon <- SplineFun(Dates)
SplineFitSAMRecon[43300:43829] = NA 

# fill missing values first
res_data_detrended$ts = 1:nrow(res_data)
res_data_detrended$cs = 1
# SOI.nc is perfectly colinear to SOInew, take it out
res_data_detrended = res_data_detrended[, names(res_data_detrended) != "SOI.nc"]
amelia_mod = amelia(x = res_data_detrended,
                    idvars = c("Month", "Year", "Date", "Doy"), # take out ts and cs from idvars
                    ts = "ts", polytime = NULL, 
                    cs = "cs",
                    parallel = "multicore", ncpus = 10,
                    #lags = names(res_data_detrended)[3:18],
                    # leads = names6,
                    m = 5)
summary(amelia_mod)
summary(amelia_mod$imputations)

plot(amelia_mod$imputations[[1]]$Date, amelia_mod$imputations[[1]]$SAMRecon)
lines(amelia_mod$imputations[[1]]$Date, amelia_mod$imputations[[1]]$SAMmonthly, col = 2)

# smooth interpolation
Dates <- seq.Date(lubridate::ymd("1900-01-01"), lubridate::ymd("2019-12-31"), by = 1)
cov_daily = data.frame(Date = Dates)

res_data_filled = amelia_mod$imputations[[1]]
# res_data_filled$SAMRecon[seq(1429, 1440, by = 3)] = 0
for (i in c(3:17)) {
  smooth = ksmooth(res_data_filled$Date, res_data_filled[,i], kernel = "normal", bandwidth = 30.4)$y
  SplineFun <- splinefun(x = res_data_filled$Date, y = smooth)
  SplineFit <- SplineFun(Dates)
  cov_daily[, names(res_data_filled)[i]] = SplineFit
}
head(cov_daily)
summary(cov_daily)

# sie_dat = cbind(data.frame(Doy, Date = sie_date),
#                 sie_detrended)
# cor_data = full_join(cov_daily, sie_dat, by = "Date")
# this.cor = cor(cor_data[,-c(1)], use = "pairwise")
# heatmap(this.cor, Rowv = NA, Colv = NA)
# round(this.cor[,17:21], 2)

# # Should we take out some of these?
# # I think AMOunsmoothed
# cov_daily = cov_daily[!(names(cov_daily) %in% c("AMOunsmoothedShort"))]

ind_cov_daily = cov_daily
# cov = cov_daily

# plot the eigenvectors
eigen_cov_ind = prcomp(ind_cov_daily[complete.cases(ind_cov_daily),-1], scale. = FALSE, center = TRUE)
melted_eigenvec <- reshape2::melt(eigen_cov_ind$rotation)
ggplot(data = melted_eigenvec, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() +
  xlab("") + ylab("") +
  scale_fill_gradient2("Loading", high = "firebrick", low = "steelblue", limits = c(-0.9, 0.9), oob = scales::squish) +
  geom_text(aes(x = Var2, y = Var1, label = round(value, 2)))

# hold out SAM and other indices for model validation?




