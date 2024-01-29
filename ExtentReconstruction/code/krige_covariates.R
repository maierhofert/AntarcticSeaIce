library("dplyr")
library("gstat")
# read all station data
station_data = read.csv("../data/fogt_predictor_data/station_data.csv")
station_data[station_data$ID == "890630_tmp", c("longname", "name", "lat", "lon")] =
  station_data[station_data$ID == "890630_SLP", c("longname", "name", "lat", "lon")]
head(station_data)
station_data = station_data[station_data$start < 1979,]

library(sp)
# extract lon lat for staton data
dsp <- SpatialPoints(station_data[,6:5], proj4string=CRS("+proj=longlat"))
dsp <- SpatialPointsDataFrame(dsp, station_data)
head(dsp)
spplot(dsp, 1)

my_stations = SpatialPoints(station_data[,6:5], proj4string=CRS("+proj=longlat"))
proj4string(my_stations)
CRS(proj4string(my_stations))
station_ID = station_data$ID

# read in weather data
all_dat = readRDS("../data/complete_dataset/all_data.RDS")
# subset to before 2021
all_dat = all_dat[all_dat$Year <= 2020,]

all_dat$tdate = all_dat$Year + all_dat$Month / 12 - 0.5
my_data = all_dat

my_dates = my_data$tdate
my_yearmon_dates = zoo::as.yearmon(my_data$tdate)
my_dates = zoo::as.Date(my_yearmon_dates)


plot(all_dat$'931190_SLP', all_dat$'948680_SLP')

used_stations = station_data$ID[station_data$ID %in% names(all_dat)]
tmp_stat = (substr(used_stations, 8, 10) == "tmp")
used_stations[tmp_stat] %>% duplicated() %>% sum()
SLP_stat = (substr(used_stations, 8, 10) == "SLP")
used_stations[SLP_stat] %>% duplicated() %>% sum()

my_stations = my_stations[station_data$ID %in% names(all_dat)]

station_match = match(used_stations, names(all_dat))

my_cov = all_dat[,station_match]
tmp = (substr(names(my_cov), 8, 10) == "tmp")
names(my_cov[,tmp]) %>% duplicated() %>% sum()
SLP = (substr(names(my_cov), 8, 10) == "SLP")
names(my_cov[,SLP]) %>% duplicated() %>% sum()

library("spacetime")
proj4string(my_stations) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
my_rural_tmp = STFDF(my_stations[tmp_stat], my_dates, 
                     data.frame(temp = as.vector(t(my_cov[,tmp]))))
summary(my_rural_tmp)
proj4string(my_rural_tmp) = CRS(proj4string(my_rural_tmp))

my_rural_SLP = STFDF(my_stations[SLP_stat], my_dates, 
                     data.frame(SLP = as.vector(t(my_cov[,SLP]))))
summary(my_rural_SLP)
proj4string(my_rural_SLP) = CRS(proj4string(my_rural_SLP))


station_ID = station_data$ID[station_data$ID %in% names(all_dat)]
# cbind(station_ID[tmp_stat],
#       my_stations[tmp_stat]@coords)


#  # subset for testing
my_rr <- my_rural_tmp[,"1898/1903"]
my_rr <- my_rural_tmp[,"1898/2020-07"]
# my_rr <- my_rural_tmp
my_rr <- as(my_rr,"STSDF")

# there should not be two stations in exactly the same location
my_rr@sp@coords %>% duplicated() %>% sum()


# rough grid
my_x1 <- seq(from=-160,to=160,by=40)
my_x2 <- seq(from=-80,to=-10,by=30)

# fine grid
my_x1 <- seq(from=-170,to=170,by=20)
my_x2 <- seq(from=-85,to=-10,by=15)

# super fine grid
my_x1 <- seq(from=-177.5,to=177.5,by=5)
my_x2 <- seq(from=-85,to=-5,by=5)


my_DE_gridded <- SpatialPoints(cbind(rep(my_x1,length(my_x2)), 
                                     rep(my_x2,each=length(my_x1))), 
                               proj4string=CRS(proj4string(my_rr@sp)))
gridded(my_DE_gridded) <- TRUE
proj4string(my_DE_gridded)

my_DE_pred <- STF(sp=as(my_DE_gridded,"SpatialPoints"), 
                  time = my_rr@time)
# time=my_rr@time[rep(c(T, F), nrow(my_rr@time) / 2),])
proj4string(my_DE_pred)

library("gstat")
# MH also compute auto-correlation function
# then compare
my_Vgm = variogramST(temp ~ 1,data = my_rr, assumeRegular = TRUE, 
                     tlags = 0:3,
                     cutoff = 2100, width = 300)
summary(my_Vgm)
plot(my_Vgm)
saveRDS(my_Vgm, "../output/my_Vgm.RDS")
# my_Vgm = readRDS("../output/my_Vgm.RDS")

# png("../plots/covariates/krige/vgm_TMP.png", height = 360, width = 540)
# plot(my_Vgm)
# dev.off()

# # there does not seem to be pronounced spatial anisotropy
# # spatial only
# my_vg = variogram(temp ~ 1,data = my_rr[,"2019-10"],
#                   alpha = c(0, 90),
#                   covariogram = TRUE)
# plot(my_vg, ylim = c(0, 2))
# 
# TheModel=vgm(psill = 0.5, model = "Mat", range = 1500,
#              kappa = 3/2, nugget = 0.2, anis = c(0, 0.5))
# FittedModel <- fit.variogram(my_vg, model=TheModel)
# ## plot results:
# plot(my_vg, model=FittedModel, as.table=TRUE)

estiStAni(my_Vgm, c(10,150))

separableModel <- vgmST("separable", 
                        # method = "Nelder-Mead", # no lower & upper needed
                        space=vgm(psill = 0.5, model = "Mat", range = 1500, 
                                  kappa = 3/2, nugget = 0.2),
                        time =vgm(psill = 0.5, model = "Mat", range = 50, 
                                  kappa = 3/2, nugget = 0.5),
                        stAni = 20,
                        sill=0.9)
# why does the anisotropy not get copied over?
my_separableVgm = fit.StVariogram(my_Vgm, separableModel,
                                  method="L-BFGS-B")
my_separableVgm
str(my_separableVgm)
plot(my_Vgm, my_separableVgm, wireframe = F, all = T)
plot(my_Vgm, my_separableVgm, wireframe = F, diff = TRUE, all = TRUE)

plot(my_separableVgm$space, cutoff = 2000)
plot(my_separableVgm$time, cutoff = 100)


# png("../plots/covariates/krige/vgm_TMP_fitted.png", height = 360, width = 720)
# plot(my_Vgm, my_separableVgm, wireframe = F, all = T)
# dev.off()

sumMetricModel <- vgmST("sumMetric",
                        space=vgm(psill = 0.5, model = "Mat", range = 1500,
                                  kappa = 3/2, nugget = 0.1, anis = c(90, 0.3)),
                        time =vgm(psill = 0.5, model = "Mat", range = 100,
                                  kappa = 3/2, nugget = 0.1),
                        joint = vgm(0.5, "Sph", 100, 0.1),
                        stAni = 20)
my_sumMetricVgm = fit.StVariogram(my_Vgm, sumMetricModel,
                                  method="L-BFGS-B")
my_sumMetricVgm
plot(my_Vgm, my_sumMetricVgm, wireframe = F, all = T)
plot(my_Vgm, my_sumMetricVgm, wireframe = F, diff = TRUE, all = TRUE)

plot(my_sumMetricVgm$space, cutoff = 2000)
plot(my_sumMetricVgm$time, cutoff = 100)
plot(my_sumMetricVgm$joint, cutoff = 10)


# somehow it is neccesary to delete this attribute in order for te krigeST function to work
attr(my_sumMetricVgm, "spatial unit") <- NULL
attr(my_separableVgm, "spatial unit") <- NULL

# do the actual kriging
my_DE_kriged <- krigeST(temp~1, data=my_rr, newdata=my_DE_pred,
                        # modelList=my_sumMetricVgm, nmax = 50
                        modelList=my_separableVgm, nmax = 50, stAni = 20)
gridded(my_DE_kriged@sp) <- TRUE
# stplot(my_DE_kriged)
# saveRDS(my_DE_kriged, "../output/my_DE_kriged_fineGrid.RDS")

saveRDS(my_DE_kriged, "../output/my_DE_kriged.RDS")
# my_DE_kriged = readRDS("../output/my_DE_kriged.RDS")



library("ggplot2")
test3 = data.frame(var = my_DE_kriged@data[, 1], 
                   lon = my_DE_kriged@sp@coords[,1],
                   lat = my_DE_kriged@sp@coords[,2],
                   date = rep(spacetime:::trimDates(my_DE_kriged), 
                              each = nrow(my_DE_kriged@sp@coords)))
ant <- ggplot2::map_data(map = "world", region = ".")
ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = var), data = test3[test3$date == "1900-01",]) +
  scale_fill_distiller("TMP", palette = "Spectral", limits = c(-2, 2), oob = scales::squish) +
  facet_wrap(~date, ncol = 12) +
  theme_bw() +
  geom_polygon(aes(x = long, y = lat, group = group),
               fill = NA, colour = "darkgrey", data = ant) +
  geom_point(aes(x = lon, y = lat), data = station_data, size = 0.5, col = "blue") +
  coord_cartesian(ylim = c(-90, 0))
ggsave("../plots/covariates/krige/TMP_1993_01.png", height = 3, width = 7)
# ggsave("../plots/covariates/krige/TMP_1993_01_fineGrid.png", height = 3, width = 7)


p1 = ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = var), data = test3) +
  scale_fill_distiller(palette = "Spectral", limits = c(-2, 2), oob = scales::squish) +
  facet_wrap(~date, ncol = 12) +
  theme_bw()
  # geom_point(aes(x = lon, y = lat), data = station_data, size = 0.001, col = "black") +
  # geom_polygon(aes(x = long, y = lat, group = group),
  #              fill = NA, colour = "darkgrey", data = ant) +
  # coord_cartesian(ylim = c(-90, 0))

# ggsave("../plots/covariates/krige.jpeg", p1, "jpeg", width = 10, height = 100,
#        limitsize = FALSE)

# ggsave("../plots/covariates/krige_fineGrid.png", p1, "png", width = 10, height = 100,
#        limitsize = FALSE)
ggsave("../plots/covariates/krige.png", p1, "png", width = 10, height = 100,
       limitsize = FALSE)

# ggplot() +
#   geom_tile(aes(x = lon, y = lat, fill = var),
#             data = test3[test3$date == "1898-08",]) +
#   scale_fill_distiller(palette = "Spectral", limits = c(-2, 2), oob = scales::squish) +
#   facet_wrap(~date, ncol = 12) +
#   theme_bw() +
#   geom_point(aes(x = lon, y = lat), data = station_data)
# head(station_data)

# ##############################################################################
# krige SLP
#  # subset for testing
my_rr_SLP <- my_rural_SLP[,"1898/1903"]
my_rr_SLP <- my_rural_SLP[,"1898/2020-07"]
my_rr_SLP <- as(my_rr_SLP,"STSDF")

# there should not be two stations in exactly the same location
my_rr_SLP@sp@coords %>% duplicated() %>% sum()

# # rough grid
# my_x1 <- seq(from=-160,to=160,by=40)
# my_x2 <- seq(from=-80,to=-10,by=30)
# 
# # fine grid
# my_x1 <- seq(from=-170,to=170,by=20)
# my_x2 <- seq(from=-85,to=-10,by=15)
# 
# my_DE_gridded <- SpatialPoints(cbind(rep(my_x1,length(my_x2)), 
#                                      rep(my_x2,each=length(my_x1))), 
#                                proj4string=CRS(proj4string(my_rr_SLP@sp)))
# gridded(my_DE_gridded) <- TRUE
# proj4string(my_DE_gridded)
# 
# my_DE_pred <- STF(sp=as(my_DE_gridded,"SpatialPoints"), 
#                   time = my_rr_SLP@time)
# proj4string(my_DE_pred)

library("gstat")
my_Vgm_SLP = variogramST(SLP ~ 1,data = my_rr_SLP, assumeRegular = TRUE, 
                         tlags = 0:3,
                         cutoff = 2100, width = 300)
summary(my_Vgm_SLP)
plot(my_Vgm_SLP)
saveRDS(my_Vgm_SLP, "../output/my_Vgm_SLP.RDS")
# my_Vgm_SLP = readRDS("../output/my_Vgm_SLP.RDS")
estiStAni(my_Vgm_SLP, c(10,150))

# there does not seem to be pronounced spatial anisotropy
# spatial only
my_vg_SLP = variogram(SLP ~ 1,data = my_rr_SLP[,"2019-12"],
                  alpha = c(0, 90))
plot(my_vg_SLP, ylim = c(0, 1))

TheModel=vgm(psill = 0.5, model = "Mat", range = 1500,
             kappa = 3/2, nugget = 0.2, anis = c(90, 1))
FittedModel <- fit.variogram(my_vg_SLP, model=TheModel)
## plot results:
plot(my_vg_SLP, model=FittedModel, as.table=TRUE)


# png("../plots/covariates/krige/vgm_SLP.png", height = 360, width = 540)
# plot(my_Vgm_SLP, wireframe = F, all = T)
# dev.off()

separableModel_SLP <- vgmST("separable", 
                            # method = "Nelder-Mead", # no lower & upper needed
                            space=vgm(psill = 0.5, model = "Mat", range = 1500, 
                                      kappa = 3/2, nugget = 0.2, anis = c(90, 0.3)),
                            time =vgm(psill = 0.5, model = "Mat", range = 50, 
                                      kappa = 3/2, nugget = 0.5),
                            stAni = 35,
                            sill=0.9)
my_separableVgm_SLP = fit.StVariogram(my_Vgm_SLP, separableModel_SLP,
                                      method="L-BFGS-B")
my_separableVgm_SLP
plot(my_Vgm_SLP, my_separableVgm_SLP, wireframe = F, all = T)
plot(my_Vgm_SLP, my_separableVgm_SLP, wireframe = F, diff = TRUE, all = TRUE)

plot(my_separableVgm_SLP$space, cutoff = 2000)
plot(my_separableVgm_SLP$time, cutoff = 100)

# png("../plots/covariates/krige/vgm_SLP_fitted.png", height = 360, width = 720)
# plot(my_Vgm_SLP, my_separableVgm_SLP, wireframe = F, all = T)
# dev.off()


sumMetricModel_SLP <- vgmST("sumMetric",
                            space=vgm(psill = 0.5, model = "Mat", range = 1500,
                                      kappa = 3/2, nugget = 0.1, anis = c(90, 0.3)),
                            time =vgm(psill = 0.5, model = "Mat", range = 100,
                                      kappa = 3/2, nugget = 0.1),
                            joint = vgm(0.5, "Sph", 100, 0.1),
                            stAni = 35)
my_sumMetricVgm_SLP = fit.StVariogram(my_Vgm_SLP, sumMetricModel_SLP,
                                      method="L-BFGS-B")
my_sumMetricVgm_SLP
plot(my_Vgm_SLP, my_sumMetricVgm_SLP, wireframe = F, all = T)
plot(my_Vgm_SLP, my_sumMetricVgm_SLP, wireframe = F, diff = TRUE, all = TRUE)

plot(my_sumMetricVgm_SLP$space, cutoff = 2000)
plot(my_sumMetricVgm_SLP$time, cutoff = 100)
plot(my_sumMetricVgm_SLP$joint, cutoff = 100)


# somehow it is neccesary to delete this attribute in order for te krigeST function to work
attr(my_sumMetricVgm_SLP, "spatial unit") <- NULL
attr(my_separableVgm_SLP, "spatial unit") <- NULL

# do the actual kriging
my_DE_kriged_SLP <- krigeST(SLP~1, data=my_rr_SLP, newdata=my_DE_pred,
                            # modelList=my_sumMetricVgm, nmax = 50
                            modelList=my_separableVgm_SLP, nmax = 50, stAni = 20)
gridded(my_DE_kriged_SLP@sp) <- TRUE
stplot(my_DE_kriged_SLP)
# saveRDS(my_DE_kriged_SLP, "../output/my_DE_kriged_SLP_fineGrid.RDS")
saveRDS(my_DE_kriged_SLP, "../output/my_DE_kriged_SLP.RDS")
# my_DE_kriged_SLP = readRDS("../output/my_DE_kriged_SLP.RDS")



library("ggplot2")
test4 = data.frame(var = my_DE_kriged_SLP@data[, 1], 
                   lon = my_DE_kriged_SLP@sp@coords[,1],
                   lat = my_DE_kriged_SLP@sp@coords[,2],
                   date = rep(spacetime:::trimDates(my_DE_kriged_SLP), 
                              each = nrow(my_DE_kriged_SLP@sp@coords)))
# get world map overlay
ant <- ggplot2::map_data(map = "world", region = ".")
ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = var), data = test4[test4$date == "1993-01",]) +
  scale_fill_distiller("SLP", palette = "Spectral", limits = c(-2, 2), oob = scales::squish) +
  facet_wrap(~date, ncol = 12) +
  theme_bw() +
  geom_point(aes(x = lon, y = lat), data = station_data, size = 0.5, col = "blue") +
  geom_polygon(aes(x = long, y = lat, group = group),
               fill = NA, colour = "darkgrey", data = ant) +
  coord_cartesian(ylim = c(-90, 0))
# ggsave("../plots/covariates/krige/SLP_1993_01_fineGrid.png", height = 3, width = 7)
ggsave("../plots/covariates/krige/SLP_1993_01.png", height = 3, width = 7)

p2 = ggplot() +
  geom_tile(aes(x = lon, y = lat, fill = var), data = test4) +
  scale_fill_distiller(palette = "Spectral", limits = c(-2, 2), oob = scales::squish) +
  facet_wrap(~date, ncol = 12) +
  theme_bw()
  # geom_point(aes(x = lon, y = lat), data = station_data, size = 0.001, col = "black") +
  # geom_polygon(aes(x = long, y = lat, group = group),
  #              fill = NA, colour = "darkgrey", data = ant) +
  # coord_cartesian(ylim = c(-90, 0))
# ggsave("../plots/covariates/krige_SLP_fineGrid.png", p2, "png", width = 10, height = 100,
#        limitsize = FALSE)
ggsave("../plots/covariates/krige_SLP.png", p2, "png", width = 10, height = 100,
       limitsize = FALSE)

