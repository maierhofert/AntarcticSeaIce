# read all data
files_mon_full = list.files("../data/fogt_predictor_data/monthly_predictor_data/short_term_all/",
                            full.names = TRUE, pattern = ".nc")
files_mon = list.files("../data/fogt_predictor_data/monthly_predictor_data/short_term_all/",
                       full.names = FALSE, pattern = ".nc")
files_mon_full = files_mon_full
files_mon = files_mon

files_mon_stations = strsplit(files_mon, c(".nc")) %>% unlist()

data_list_mon = list()

# create data frame for station data
station_data = data.frame(ID = files_mon_stations)
for (i in 1:length(files_mon)) {
  this.file = files_mon_full[i]
  this.nc = nc_open(this.file)
  station_data$longname[i] = this.nc$var[[1]]$longname
}
station_data
write.csv2(station_data, "../data/fogt_predictor_data/station_data_raw.csv")

# manually fill in station name and location from longname
station_data = read.csv("../data/fogt_predictor_data/station_data.csv")
station_data$var = sapply(station_data$ID, function(x) {
  strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])]
})

# ##############################################################################
library("ggplot2")
# plot map of antarctica
# taken from ANNCYCLEPlotFutureIceOutlines.R
ant <- ggplot2::map_data(map = "world", region = ".",
                         # orientation = c(-90, 0, 0),
                         wrap = c(-180, 180, -90))

# plot a map of antarctica
antmap <- ggplot(ant, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white", colour = "black") +
  coord_map("ortho", orientation=c(-90, 0, 0), ylim = c(-90,-50))
# Plot it in cartesian coordinates
antmap

# geographic projection
density_dat = density(station_data[!duplicated(station_data$name) &
                                     station_data$start < 1980,"lon"],
                      bw = 10,
                      from = -180, to = 180)
density_dat = data.frame(x = density_dat$x, y = density_dat$y)

density_dat_lat = density(station_data[!duplicated(station_data$name) &
                                     station_data$start < 1980,"lat"],
                      bw = 5,
                      from = -90, to = 0)
density_dat_lat = data.frame(x = density_dat_lat$x, y = density_dat_lat$y)
density_dat_lat = density_dat_lat[order(density_dat_lat$x),]
plot(density_dat_lat$y, density_dat_lat$x, type = "l")

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
  geom_point(aes(x = lon, y = lat),
             data = station_data[!duplicated(station_data$name) &
                                   station_data$start < 1980,]) +
  ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name),
                           max.overlaps = 20,
                           data = station_data[!duplicated(station_data$name) &
                                                 station_data$start < 1980,]) +
  geom_line(aes(x = x, y = y*2000 - 85), data = density_dat) +
  geom_path(aes(x = y*2000 - 180, y = x), data = density_dat_lat)
ggsave(filename = paste0("../plots/Antarctica/AntarcticaLonLat_stationMonthly.png"),
       width = 14, height = 7,
       device = png())
dev.off()


################################################################################
# daily stations
station_data_daily = readxl::read_xlsx("../data/fogt_predictor_data/daily_predictor_data/daily_GHCND_data/Key.xlsx",
                                       col_names = c("ID", "lat", "lon", "var", "start", "end", "name", "ID2"))

head(station_data_daily)

all_stations = full_join(station_data[-1], station_data_daily,
                         by = c("ID", "lon", "lat", "name", "var", "start", "end"))
saveRDS(all_stations, "../data/fogt_predictor_data/all_stations.RDS")

# ##############################################################################
# plot map of antarctica
# taken from ANNCYCLEPlotFutureIceOutlines.R

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
  geom_point(aes(x = lon, y = lat),
             data = station_data_daily[!duplicated(station_data_daily$name) &
                                         station_data_daily$start < 1980,]) +
  ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name),
                           max.overlaps = 20,
                           data = station_data_daily[!duplicated(station_data_daily$name) &
                                                       station_data_daily$start < 1980,])
ggsave(filename = paste0("../plots/Antarctica/AntarcticaLonLat_Stationdaily.png"),
       width = 14, height = 7,
       device = png())
dev.off()

