# combine station location and variable importance from RF
mod_list = readRDS("../data/models/mod_list_RF_all.RDS")
res_data = readRDS("../data/models/res_data_RF_all.RDS")

all_stations = readRDS("../data/fogt_predictor_data/all_stations.RDS")

# check out variable importance
png(paste0("../plots/performance/VarImp_allCovLag_", res_data$sector[1], ".png"),
    height = 480, width = 480)
varImpPlot(mod_list[[1]], main = res_data$pretty_sector[1])
dev.off()
png(paste0("../plots/performance/VarImp_allCovLag_", res_data$sector[2], ".png"),
    height = 480, width = 480)
varImpPlot(mod_list[[2]], main = res_data$pretty_sector[2])
dev.off()
png(paste0("../plots/performance/VarImp_allCovLag_", res_data$sector[3], ".png"),
    height = 480, width = 480)
varImpPlot(mod_list[[3]], main = res_data$pretty_sector[3])
dev.off()
png(paste0("../plots/performance/VarImp_allCovLag_", res_data$sector[4], ".png"),
    height = 480, width = 480)
varImpPlot(mod_list[[4]], main = res_data$pretty_sector[4])
dev.off()
png(paste0("../plots/performance/VarImp_allCovLag_", res_data$sector[5], ".png"),
    height = 480, width = 480)
varImpPlot(mod_list[[5]], main = res_data$pretty_sector[5])
dev.off()

(importance(mod_list[[1]]) %>% 
    rownames() %>%
    substring(first = 2)) ==
  c("random", names(cov))

dim(cov)

RFi_importance = function(i) {
  this.importance = importance(mod_list[[i]])
  ret = data.frame(sector = res_data$sector[i],
                   pretty_sector = res_data$pretty_sector[i],
                   year_block_ind = res_data$year_block_ind[i],
                   ID = rownames(this.importance) %>%
                     substring(first = 2), 
                   IncNodePurity = this.importance)
  ret$ID[2] = "Doy"
  ret
}
RF_importance = RFi_importance(1)
for (i in 2:nrow(res_data)) {
  RF_importance = rbind(RF_importance, RFi_importance(i))
}

RF_importance$match_ID = sapply(RF_importance$ID, function(x) {
  strsplit(x, split = "_l")[[1]][1]
})
# RF_importance$match_ID = sapply(RF_importance$match_ID, function(x) {
#   x = strsplit(x, split = c("TMIN_"))[[1]][length(strsplit(x, split = c("TMIN_"))[[1]])]
#   x = strsplit(x, split = c("TMAX_"))[[1]][length(strsplit(x, split = c("TMAX_"))[[1]])]
#   x
# })

RF_importance$lag = sapply(RF_importance$ID, function(x) {
  strsplit(x, split = "_l")[[1]][2]
})
RF_importance$lag[is.na(RF_importance$lag)] = 0

# combine with station data
all_stations$match_ID = all_stations$ID
daily_ind = sapply(all_stations$ID, function(x) {
  length(strsplit(x, "_")[[1]]) == 1
})
all_stations$match_ID[daily_ind] = paste(all_stations$var, all_stations$ID, 
                                                        sep = "_")[daily_ind]
# add a fake station for day of year
doy_station = data.frame(ID = "00", longname = "day of year", name = "day of year", 
                lat = -85, lon = 120, ID.1 = "", 
                start = "1900", end = "2019", var = "tmp", ID2 = "",
                match_ID = "Doy")
all_stations = rbind(doy_station, all_stations)
# combined data fo RF importance and stations
station_importance = left_join(RF_importance, all_stations, by = c("match_ID" = "match_ID"))
table(station_importance$var, useNA = "always")

# 
# ##############################################################################
# plot importance of every station
station_importance_byName = station_importance %>%
  group_by(sector, pretty_sector, year_block_ind, name, lat, lon) %>%
  summarise(IncNodePurity = sum(IncNodePurity)) %>%
  group_by(sector, year_block_ind) %>%
  filter(IncNodePurity > quantile(IncNodePurity, 0.75)) %>%
  ungroup()

# plot map of antarctica
# taken from ANNCYCLEPlotFutureIceOutlines.R
ant <- ggplot2::map_data(map = "world", region = ".", 
                         # orientation = c(-90, 0, 0),
                         wrap = c(-180, 180, -90))

filter_dat = filter(station_importance_byName, 
                    sector == "King_Hakon",
                    year_block_ind == 1)
# geographic projection
ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group),
               data = ant,
               fill = "gray70", colour = "gray70") +
  # coord_map("stereographic", orientation=c(-90, 0, 0), ylim = c(-90,-50)) +
  scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
  scale_y_continuous(name = "latitude") +
  coord_cartesian(ylim = c(-90, 0)) +
  # sea ice based sectors
  geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
  # sea ice based sectors
  geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5), y = -3,
                               label = c("East\nAntarctica",
                                         "Ross Sea",
                                         "Amundsen\nBellingsh.\nSea",
                                         "Weddell\nSea",
                                         "King Hakon")),
             aes(x = x, y = y, label = label),
             color = "steelblue") +
  theme_minimal() +
  # add stations
  geom_point(aes(x = lon, y = lat, size = IncNodePurity, alpha = IncNodePurity), 
             data = filter_dat) +
  ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = IncNodePurity), 
                           max.overlaps = 20,
                           data = filter_dat) +
  ggtitle(paste0(filter_dat$pretty_sector[1], ", CV-fold: ", filter_dat$year_block_ind[1], ", 25% most important stations")) +
  scale_alpha("Variable\nImportance") +
  scale_size("Variable\nImportance")

ggsave(filename = paste0("../plots/Antarctica/AntarcticaLonLat_stationImp_",
                         filter_dat$sector[1], "_", filter_dat$year_block_ind[1], ".png"),
       width = 14, height = 7,
       device = png())
dev.off()


# ###############################################################################
# below is just for plotting

# geographic projection
base_plot = ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group),
               data = ant,
               fill = "gray70", colour = "gray70") +
  # coord_map("stereographic", orientation=c(-90, 0, 0), ylim = c(-90,-50)) +
  scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
  scale_y_continuous(name = "latitude") +
  coord_cartesian(ylim = c(-90, 0)) +
  # sea ice based sectors
  geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
  # sea ice based sectors
  geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5), y = -3,
                               label = c("East\nAntarctica",
                                         "Ross Sea",
                                         "Amundsen\nBellingsh.\nSea",
                                         "Weddell\nSea",
                                         "King Hakon")),
             aes(x = x, y = y, label = label),
             color = "steelblue") +
  theme_minimal()
for (i in which(res_data$year_block_ind %in% c(1, 8))) {
  filter_dat = filter(station_importance_byName, 
                      sector == res_data$sector[i],
                      year_block_ind == res_data$year_block_ind[i])
  
  # add stations
  base_plot +
    geom_point(aes(x = lon, y = lat, size = IncNodePurity, alpha = IncNodePurity), 
               data = filter_dat) +
    ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = IncNodePurity), 
                             max.overlaps = 20,
                             data = filter_dat) +
    ggtitle(paste0(filter_dat$pretty_sector[1], ", CV-fold: ", filter_dat$year_block_ind[1], ", 25% most important stations")) +
    scale_alpha("Variable\nImportance") +
    scale_size("Variable\nImportance")
  
  ggsave(filename = paste0("../plots/Antarctica/RF_StationImp/AntarcticaLonLat_stationImp_",
                           filter_dat$sector[1], "_", filter_dat$year_block_ind[1], ".png"),
         width = 14, height = 7,
         device = png())
  dev.off()
}

# ##############################################################################
# split it up by lag
station_importance_byNameLag = station_importance %>%
  group_by(sector, pretty_sector, year_block_ind, name, lat, lon, lag) %>%
  summarise(IncNodePurity = sum(IncNodePurity)) %>%
  group_by(sector, year_block_ind) %>%
  filter(IncNodePurity > quantile(IncNodePurity, 0.90)) %>%
  ungroup()

for (i in which(res_data$year_block_ind %in% c(1, 8))) {
  
  filter_dat = filter(station_importance_byNameLag, 
                      sector == res_data$sector[i],
                      year_block_ind == res_data$year_block_ind[i])
  
  this.base = ggplot(filter_dat) +
    geom_polygon(aes(x = long, y = lat, group = group),
                 data = ant,
                 fill = "gray70", colour = "gray70") +
    # coord_map("stereographic", orientation=c(-90, 0, 0), ylim = c(-90,-50)) +
    scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(name = "latitude") +
    coord_cartesian(ylim = c(-90, 0)) +
    # sea ice based sectors
    geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
    # sea ice based sectors
    geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5), y = -3,
                                 label = c("East\nAntarctica",
                                           "Ross Sea",
                                           "Amundsen\nBellingsh.\nSea",
                                           "Weddell\nSea",
                                           "King Hakon")),
               aes(x = x, y = y, label = label),
               color = "steelblue") +
    theme_minimal()
  
  # add stations
  this.base +
    facet_wrap(vars(factor(as.character(lag), levels = c("0", "7", "29", "61")))) +
    geom_point(aes(x = lon, y = lat, size = IncNodePurity, alpha = IncNodePurity)) +
    ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = IncNodePurity), 
                             max.overlaps = 20) +
    ggtitle(paste0(filter_dat$pretty_sector[1], ", CV-fold: ", filter_dat$year_block_ind[1], ", 10% most important stations")) +
    scale_alpha("Variable\nImportance") +
    scale_size("Variable\nImportance")
  
  ggsave(filename = paste0("../plots/Antarctica/RF_StationImp_lag/AntarcticaLonLat_stationImp_",
                           filter_dat$sector[1], "_", filter_dat$year_block_ind[1], "_lag.png"),
         width = 14, height = 7,
         device = png())
  dev.off()
}

# ###########################################################################
# ##############################################################################
# split it up by type of variable
table(station_importance$var, useNA = "always")
station_importance_byNameVar = station_importance %>%
  group_by(sector, pretty_sector, year_block_ind, name, lat, lon, var) %>%
  summarise(IncNodePurity = sum(IncNodePurity)) %>%
  group_by(sector, year_block_ind) %>%
  filter(IncNodePurity > quantile(IncNodePurity, 0.90)) %>%
  ungroup()

for (i in which(res_data$year_block_ind %in% c(1, 8))) {
  
  filter_dat = filter(station_importance_byNameVar, 
                      sector == res_data$sector[i],
                      year_block_ind == res_data$year_block_ind[i])
  
  this.base = ggplot(filter_dat) +
    geom_polygon(aes(x = long, y = lat, group = group),
                 data = ant,
                 fill = "gray70", colour = "gray70") +
    # coord_map("stereographic", orientation=c(-90, 0, 0), ylim = c(-90,-50)) +
    scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(name = "latitude") +
    coord_cartesian(ylim = c(-90, 0)) +
    # sea ice based sectors
    geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
    # sea ice based sectors
    geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5), y = -3,
                                 label = c("East\nAntarctica",
                                           "Ross Sea",
                                           "Amundsen\nBellingsh.\nSea",
                                           "Weddell\nSea",
                                           "King Hakon")),
               aes(x = x, y = y, label = label),
               color = "steelblue") +
    theme_minimal()
  
  # add stations
  this.base +
    facet_wrap(vars(var)) +
    geom_point(aes(x = lon, y = lat, size = IncNodePurity, alpha = IncNodePurity)) +
    ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = IncNodePurity), 
                             max.overlaps = 20) +
    ggtitle(paste0(filter_dat$pretty_sector[1], ", CV-fold: ", filter_dat$year_block_ind[1], ", 10% most important stations")) +
    scale_alpha("Variable\nImportance") +
    scale_size("Variable\nImportance")
  
  ggsave(filename = paste0("../plots/Antarctica/RF_StationImp_var/AntarcticaLonLat_stationImp_",
                           filter_dat$sector[1], "_", filter_dat$year_block_ind[1], "_var.png"),
         width = 14, height = 7,
         device = png())
  dev.off()
}








