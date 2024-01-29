# this file reconstructs the analysis done by Fogt et al 2022
# this time we do the cross validation over the entire variable selection process
library("dplyr")
library("ggplot2")

# read all_data_s created in fogt_rebuild.R
all_data_s = readRDS("all_data_seasonal_lagged.RDS")
all_dat_s = all_data_s[between(all_data_s$season_year, 1979, 2020),]

# helper function, coefficient of efficiency
my_CE = function(truth, pred) {
  1 - sum((truth - pred)^2) / sum((truth - mean(truth))^2)
}

# ##############################################################################
# set up a nice cross validation
# sectors
sector = colnames(all_data_s)[3:8]
pretty_sector = c("Total", "King Hakon VII", "Ross Sea", "East Antarctica", 
                  "Weddell Sea", "Bellingshausen\nAmundsen Sea")
season = c("DJF", "MAM", "JJA", "SON")
n_season = 1:4
cor_threshold = c(0.1, 0.05, 0.025, 0.01)
n_indices = 0:9
n_pca = seq(2, 35, by = 2)


# set up results data
res_data = expand.grid(sector = sector,
                       n_season = n_season,
                       cor_threshold = cor_threshold,
                       n_indices = n_indices,
                       n_pca = n_pca)
# add all observations regardless of significance
res_data_allObs = expand.grid(sector = sector,
                              n_season = n_season,
                              cor_threshold = 1,
                              n_indices = 9,
                              n_pca = n_pca)
res_data = rbind(res_data, res_data_allObs)

# add sector 
res_data$sector = as.character(res_data$sector)
res_data$pretty_sector = factor(pretty_sector, levels = pretty_sector)
# add season
res_data$season = "DJF"
res_data$season[res_data$n_season == 2] = "MAM"
res_data$season[res_data$n_season == 3] = "JJA"
res_data$season[res_data$n_season == 4] = "SON"
res_data$season = factor(res_data$season, 
                         levels = c("DJF", "MAM", "JJA", "SON"))

# performance measures
res_data$val_cor = NA
res_data$val_ce = NA
res_data$val_rmse = NA
res_data$rmse0 = NA

# number of models per sector and season
# divide by number of retained PC CV
sum(res_data$sector == "total" & res_data$n_season == 1) / length(n_pca)

# we are going through 16728 models
nrow(res_data)
# it should be 28.8k models
5*4*4*10*36
sector = colnames(all_data_s)[3:8]
season = c("DJF", "MAM", "JJA", "SON")
n_season = 1:4
cor_threshold = c(0.1, 0.05, 0.025, 0.01)
n_indices = 0:9
n_pca = seq(2, 35, by = 2)

# add covariate prediction strength
all_cov_names = names(all_dat_s)[9:664]
res_data$intercept = 0
res_data[,all_cov_names] = 0

for (j in 1:nrow(res_data)) {
  if (j %% 100 == 0) print(j)
  # what to use for reconstruction
  # input variables
  # which sector
  this.sector = res_data$sector[j]
  # which season
  this.season = res_data$season[j]
  this.n_season = res_data$n_season[j]
  # which correlation significance threshold to use   
  this.cor_threshold = res_data$cor_threshold[j]
  # number of climate indices to include
  this.n_indices = res_data$n_indices[j]
  # how many PC should be retained
  this.n_pca = res_data$n_pca[j]
  
  yrs = 1979:2020
  preds = rep(NA, length(yrs))
  my_pred = rep(NA, length(yrs))
  # read out the coefficients on original level
  coefs = matrix(0, nrow = length(yrs), ncol = 1 + length(all_cov_names))
  colnames(coefs) = c("intercept", all_cov_names)
  for (i in 1:length(yrs)) {
    this.pred_year = yrs[i]
    # delete two years before and after the year to be predicted
    del_years = -2:2 + this.pred_year
    this.all_dat_s = all_dat_s[!all_dat_s$season_year %in% del_years,]
    
    # computed variables
    # compute correlation
    # compute separate correlation matrices per season
    this.corr_test = psych::corr.test(this.all_dat_s[this.all_dat_s$n_season == this.n_season, this.sector],
                                      this.all_dat_s[this.all_dat_s$n_season == this.n_season, 9:655],
                                      use = "pairwise.complete.obs")
    this.corr_p = this.corr_test$p
    
    # which covariates meet the correlation significance threshold in this season
    this.covariates = colnames(this.corr_p)[this.corr_p < this.cor_threshold]
    # sort out indices that should not be used
    indices = c("IPOunfilteredV5",
                "AMOunsmoothedLong",
                "PDOnew",
                "SOInew",
                "SSTsNino1.2",
                "SSTsNino3.4",
                "SSTsNino3Mean",
                "SSTsNino4Mean",
                "SAMRecon")
    if (this.n_indices < 9) {
      # delete unwanted climate indices from covariates
      this.del_indices = indices[(this.n_indices + 1):9]
      for (ind in this.del_indices) {
        del = grep(ind, this.covariates)
        if (length(del) > 0) {this.covariates = this.covariates[-del]}
      }
    }
    
    if (length(this.covariates) >= 2) {
      # do a PC transform
      this.pca = prcomp(this.all_dat_s[this.all_dat_s$season == this.season, this.covariates], center = FALSE)
      # which PCs are the most strongly correlated?
      this.pc_cor = cor(this.all_dat_s[this.all_dat_s$season == this.season, this.sector],
                        this.pca$x[,1:min(length(this.covariates), 35)])
      # order the PC by correlation magnitude
      this.strongest_pc = colnames(this.pc_cor)[order(abs(this.pc_cor), decreasing = TRUE)]
      # get the strongest corelations
      this.strongest_pc = this.strongest_pc[1:min(this.n_pca, length(this.covariates))]
      
      # only retain first 35 principal components
      this.mod_data = data.frame(this.all_dat_s[this.all_dat_s$season == this.season, 
                                                1:8],
                                 this.pca$x[,this.strongest_pc])
      this.newdata = predict(this.pca, 
                             all_dat_s[all_dat_s$season_year == this.pred_year &
                                         all_dat_s$season == this.season,]) %>%
        data.frame()
      # create formula
      this.formula = paste0(this.sector,
                            " ~ ",
                            paste(this.strongest_pc,
                                  collapse = " + ", sep = "")) %>%
        formula()
    } else if (length(this.covariates) == 1) {
      this.mod_data = this.all_dat_s[this.all_dat_s$season == this.season,]
      this.newdata = all_dat_s[all_dat_s$season_year == this.pred_year & all_dat_s$season == this.season,]
      # create formula
      this.formula = paste0(this.sector,
                            " ~ ",
                            paste(this.covariates,
                                  collapse = " + ", sep = "")) %>%
        formula()
    } else {
      stop("This is not implemented")
    }
    
    this.mod = lm(this.formula, 
                  data = this.mod_data)
    # mod_list[[i]] = this.mod
    # summary(this.mod)
    preds[i] = predict(this.mod, 
                       newdata = this.newdata)
    # this.coef %*% c(1, unlist(this.newdata[,names(this.coef[-1])]))
    
    # compute regression coefficients for original data
    # intercept stays the same
    this.coef = coefficients(this.mod)
    # multiply PC coefficients by rotation to get implicit original coefficients
    new_coef = this.pca$rotation[,this.strongest_pc] %*% this.coef[-1]
    # this is the data for the left out year
    newdata = all_dat_s[all_dat_s$season_year == this.pred_year &
                          all_dat_s$season == this.season, this.covariates] %>% as.matrix(nrow = 1)
    # this is exactly the same prediction, but using the original data as opposed to the PCs
    my_pred[i] = this.coef[1] + newdata %*% new_coef
    # add coefficients for this year to our model
    coefs[i,"intercept"] = this.coef[1]
    coefs[i,this.covariates] = new_coef
  }
  # compute average strength of covariates and add it to the results
  res_data[j, c("intercept", all_cov_names)] = colMeans(abs(coefs))
  # # plot predicted vs fitted
  # plot(1979:2020, this.mod_data[,this.sector], type = "b", col = "green")
  # abline(h = 0, lty = 2)
  # points(1979:2020, preds, type = "b")
  
  res_data[j,"val_cor"] = cor(preds, all_dat_s[all_dat_s$season == this.season, this.sector])
  res_data[j,"val_ce"] = my_CE(all_dat_s[all_dat_s$season == this.season, this.sector], preds)
  res_data[j,"val_rmse"] = sqrt(mean((preds - all_dat_s[all_dat_s$season == this.season, this.sector])^2))
  res_data[j,"rmse0"] = sqrt(mean((all_dat_s[all_dat_s$season == this.season, this.sector])^2))
}


# save version with duplicates and n_pca
res_data_old = res_data

# find the optimal number of principal components to be retained
res_data3 = res_data %>%
  group_by(sector, n_season, cor_threshold, n_indices) %>%
  slice(which.max(val_ce))
table(res_data3$n_pca)
duplicated(res_data3[, c("sector", "n_season", "val_cor", "val_ce", "val_rmse")]) %>% sum()

# remove duplicated models
res_data4 = res_data3[!duplicated(res_data3[, c("sector", "n_season", "val_cor", "val_ce", "val_rmse")]),]
table(res_data4$n_pca)
hist(res_data4$n_pca)

# reduced data to only unique models
res_data_red = res_data4
hist(res_data_red$val_ce[res_data_red$val_ce > 0])
hist(res_data_red$val_cor)

# save without duplicates
# saveRDS(res_data_red, "fogt_rebuild_results_fixLOOCV.RDS")
# res_data_red = readRDS("../output/fogt_rebuild_results_fixLOOCV.RDS")
saveRDS(res_data_red, "fogt_rebuild_results_fixLOOCV_lindetrend.RDS")
# res_data_red = readRDS("fogt_rebuild_results_fixLOOCV_lindetrend.RDS")


# combine results with results from original method
res_data_red$method = "complete LOOCV"
# read in original methods results
# res_data_red_og = readRDS("fogt_rebuild_results.RDS")
res_data_red_og = readRDS("fogt_rebuild_results_lindetrend.RDS")
res_data_red_og$method = "original"
# combine results
res_data_red = rbind(res_data_red, res_data_red_og)
res_data_red$method = factor(res_data_red$method, 
                             levels = c("original", "complete LOOCV"))
library("ggplot2")

# plot combined results
ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~pretty_sector, nrow = 1) +
  geom_boxplot(aes(x = season, y = val_cor, col = method)) +
  scale_color_manual(values = c("firebrick", "steelblue"), 
                     labels = c("original", "complete\nLOOCV")) +
  scale_x_discrete("Season") + 
  scale_y_continuous("Validation Correlation") +
  theme_bw() +
  coord_cartesian(ylim = c(-0.5, 1))
ggsave("plots/fogt_rebuild_valCor_fixLOOCV.png", height = 4, width = 12)

ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~pretty_sector, nrow = 1) +
  geom_boxplot(aes(x = season, y = val_ce, col = method)) +
  scale_color_manual(values = c("firebrick", "steelblue"),
                     labels = c("original", "complete\nLOOCV")) +
  scale_x_discrete("Sector") +
  scale_y_continuous("Validation CE") +
  theme_bw() +
  coord_cartesian(ylim = c(-0.5, 1))
ggsave("plots/fogt_rebuild_valCe_fixLOOCV.png", height = 4, width = 12)

# how many models do we have per sector and season
res_data_red %>%
  filter(sector != "Total", method == "complete LOOCV") %>%
  group_by(season, sector, method) %>%
  summarise(n())
tab = table(res_data_red[,c("sector", "season", "method")])
tab

# print out some more summary statistics
res_data_red %>%
  filter(sector != "Total", method == "complete LOOCV") %>%
  group_by(season, sector) %>%
  summarise(val_ce = mean(val_ce))

res_data_red %>%
  filter(sector != "Total", method == "complete LOOCV") %>%
  group_by(sector) %>%
  summarise(val_ce = mean(val_ce))

res_data_red %>%
  filter(sector != "Total", method == "complete LOOCV") %>%
  group_by(sector) %>%
  summarise(val_cor = mean(val_cor))

res_data_red %>%
  filter(sector != "Total", method == "complete LOOCV") %>%
  group_by(sector) %>%
  summarise(val_cor = mean(val_cor)) %>%
  summarise(val_cor = mean(val_cor))
# 0.370

res_data_red %>%
  filter(sector != "Total", method == "complete LOOCV") %>%
  group_by(sector) %>%
  summarise(val_ce = mean(val_ce)) %>%
  summarise(val_ce = mean(val_ce))
# 0.0274




# ##############################################################################
# can we look into the distribution of covariates
table(res_data_red$method)

# ##############################################################################
res_data_red = res_data_red[res_data_red$method == "complete LOOCV",]
# compute some kind of variable importance
var_imp = res_data_red[,all_cov_names]
# recombine all lags with original variable
var_imp_total = var_imp[,1:164]
# match all lags of the same variable back together
match_ID = sapply(colnames(var_imp), function(x) {
  strsplit(x, split = "_l")[[1]][1]
})
for (i in 1:164) {
  var_imp_total[,i] = rowSums(var_imp[,match_ID == match_ID[i]])
}

# read in station data
all_stations = read.csv("station_data.csv")
all_stations$matchID = paste0("V", all_stations$ID)
# just a quick check
# only indices do not have. station associated
which(!(colnames(var_imp_total) %in% all_stations$matchID))
colnames(var_imp_total)[!(colnames(var_imp_total) %in% all_stations$matchID)]

# subset to stations that are in the model
all_stations = all_stations[all_stations$matchID %in% colnames(var_imp_total),]

all_stations$var = sapply(all_stations$ID, function(x) {
  strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])]
})
# add fake stations for climate indices
index_stations = data.frame(ID = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                   "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                   "SSTsNino4Mean.1854"),
                            longname = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                         "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                         "SSTsNino4Mean.1854"),
                            name = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                     "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                     "SSTsNino4Mean.1854"),
                            lat = -85, lon = seq(-150,150, length.out = 9),
                            ID.1 = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                     "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                     "SSTsNino4Mean.1854"),
                            matchID = c("VAMOunsmoothedLong", "VIPOunfilteredV5", "VPDOnew", "VSAMRecon",
                                        "VSOInew", "VSSTsNino1.2Mean.1854", "VSSTsNino3.4Mean.1854", "VSSTsNino3Mean.1854",
                                        "VSSTsNino4Mean.1854"),
                            start = "1900", end = "2019", var = "ind")
all_stations = rbind(all_stations, index_stations)

# check that all names are alphabetically ordered
order(colnames(var_imp_total))
order(all_stations$matchID)

# add importance
station_importance = all_stations
# station_importance[,c("importance_King_Haakon", "importance_Ross",
#                       "importance_East_Antarctica", "importance_Weddell",
#                       "importance_Bellingshausen_Amundsen")] = t(var_imp_total)
station_importance[,paste0("M",1:nrow(var_imp_total))] = t(var_imp_total)

# aggregate by station
station_importance_agg = station_importance %>%
  group_by(name,  lat, lon) %>%
  summarise(across(M1:M343, sum))

# plot station importance
ant <- ggplot2::map_data(map = "world", region = ".",
                         # orientation = c(-90, 0, 0),
                         wrap = c(-180, 180, -90))
sectors = c("King_Hakon", "Ross",
            "East_Antarctica", "Weddell",
            "Bellingshausen_Amundsen_Sea")
pretty_sectors = c("King Haakon VII", "Ross Sea",
                   "East Antarctica", "Weddell Sea",
                   "Bellingshausen Amundsen Sea")
# first empty out plots
file.remove(list.files("plots/station_importance_fixLOOCV", full.names = TRUE))
for (i in 1:nrow(res_data_red)) {
  this.sector = res_data_red$sector[i]
  this.pretty_sector = res_data_red$pretty_sector[i]
  this.season = res_data_red$season[i]
  this.pretty_season = c("DJF", "MAM", "JJA", "SON")[this.season]
  # make sure we get the right importance per sector
  filter_dat = station_importance_agg[,c(1:3, 3+i)]
  names(filter_dat)[4] = "importance"
  # subset to most important 25% of stations
  filter_dat = filter_dat %>% ungroup() %>% 
    filter(importance > quantile(importance, 0.75))
  # get density without the climate indices
  station_importance_agg_noInd = station_importance_agg[!(station_importance_agg$name %in%
                                                            c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                                              "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                                              "SSTsNino4Mean.1854")),]
  density_dat = density(station_importance_agg_noInd$lon,
                        bw = 10,
                        from = -180, to = 180,
                        weights = unlist(station_importance_agg_noInd[,3+i]) / sum(unlist(station_importance_agg_noInd[,3+i])))
  # density_dat = density(filter_dat$lon,
  #                       bw = 10,
  #                       from = -180, to = 180,
  #                       weights = unlist(filter_dat[,4]) / sum(unlist(filter_dat[,4])))
  density_dat = data.frame(x = density_dat$x, y = density_dat$y)
  # geographic projection
  my_col = rep("steelblue", 5)
  my_col[which(this.sector == sectors)] = "firebrick"
  if (this.sector == "total") {
    my_col = rep("firebrick", 5)
  }
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group),
                 data = ant,
                 fill = "gray70", colour = "gray70") +
    # coord_map("stereographic", orientation=c(-90, 0, 0), ylim = c(-90,-50)) +
    scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(name = "latitude") +
    coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 0)) +
    # sea ice based sectors
    geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
    # sea ice based sectors
    geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5)[c(5,2,1,4,3)], y = -3,
                                 label = c("East\nAntarctica",
                                           "Ross Sea",
                                           "Amund.\nBel. Sea",
                                           "Weddell\nSea",
                                           "King Haakon")[c(5,2,1,4,3)]),
               aes(x = x, y = y, label = label),
               color = my_col) +
    theme_minimal() +
    # add stations
    geom_point(aes(x = lon, y = lat, size = importance, alpha = importance),
               data = filter_dat, show.legend = FALSE) +
    ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = importance),
                             max.overlaps = 20,
                             show.legend = FALSE,
                             data = filter_dat) +
    ggtitle(paste0(this.pretty_sector, ", ", this.pretty_season, ", ", "25% most important stations")) +
    scale_alpha("Variable\nImportance") +
    scale_size("Variable\nImportance") +
    geom_line(aes(x = x, y = y * 2000 - 85), data = density_dat)
  ggsave(filename = paste0("plots/station_importance_fixLOOCV/StationImp_",
                           res_data_red[i, "sector"],
                           "_ssn", res_data_red[i, "season"],
                           "_corP", res_data_red[i,"cor_threshold"],
                           "_nInd", res_data_red[i,"n_indices"],
                           "_nPC", res_data_red[i,"n_pca"], 
                           ".png"),
         width = 8, height = 5,
         device = png())
  dev.off()
  
  # filter_dat2 = slice_sample(filter_dat, n = 200, weight_by = filter_dat$importance, replace = TRUE)
  # table(filter_dat2$name)
  # 
  # ggplot() +
  #   # background density
  #   stat_density_2d(aes(x=lon, y=lat, fill = ..density..), h = c(15, 30),
  #                   data = filter_dat2, geom = "raster", contour = FALSE, show.legend = FALSE) +
  #   # scale_fill_distiller(palette="Spectral", direction=-1) +
  #   scale_fill_distiller(palette="RdBu", direction=-1) +
  #   # scale_fill_distiller(palette="Blues", direction=1) +
  #   geom_polygon(aes(x = long, y = lat, group = group),
  #                data = ant,
  #                fill = NA,
  #                colour = "gray70") +
  #   scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
  #   scale_y_continuous(name = "latitude") +
  #   coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 0)) +
  #   # sea ice based sectors
  #   geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
  #   # sea ice based sectors
  #   geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5)[c(5,2,1,4,3)], y = -3,
  #                                label = c("East\nAntarctica",
  #                                          "Ross Sea",
  #                                          "Amund.\nBel. Sea",
  #                                          "Weddell\nSea",
  #                                          "King Haakon")[c(5,2,1,4,3)]),
  #              aes(x = x, y = y, label = label),
  #              color = my_col) +
  #   theme_minimal() +
  #   # add stations
  #   geom_point(aes(x = lon, y = lat, size = importance, alpha = importance),
  #              data = filter_dat, show.legend = FALSE) +
  #   ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = importance),
  #                            max.overlaps = 20,
  #                            show.legend = FALSE,
  #                            data = filter_dat) +
  #   ggtitle(paste0(this.pretty_sector, ", ", "25% most important stations")) +
  #   scale_alpha("Variable\nImportance") +
  #   scale_size("Variable\nImportance") +
  #   geom_line(aes(x = x, y = y * 2000 - 85), data = density_dat)
  # ggsave(filename = paste0("../plots/BayesianModel/StationImpFill_",
  #                          this.sector, ".png"),
  #        width = 8, height = 5,
  #        device = png())
  # dev.off()
}

# ##############################################################################
# compute some kind of variable importance
var_imp_ss = res_data_red %>%
  group_by(sector, pretty_sector, n_season, season) %>%
  summarise(across(V619010_SLP:VSSTsNino4Mean.1854_l3, mean)) %>%
  ungroup()
var_imp = var_imp_ss[,all_cov_names]

# recombine all lags with original variable
var_imp_total = var_imp[,1:164]
# match all lags of the same variable back together
match_ID = sapply(colnames(var_imp), function(x) {
  strsplit(x, split = "_l")[[1]][1]
})
for (i in 1:164) {
  var_imp_total[,i] = rowSums(var_imp[,match_ID == match_ID[i]])
}


# read in station data
all_stations = read.csv("station_data.csv")
all_stations$matchID = paste0("V", all_stations$ID)
# just a quick check
# only indices do not have. station associated
which(!(colnames(var_imp_total) %in% all_stations$matchID))
colnames(var_imp_total)[!(colnames(var_imp_total) %in% all_stations$matchID)]

# subset to stations that are in the model
all_stations = all_stations[all_stations$matchID %in% colnames(var_imp_total),]

all_stations$var = sapply(all_stations$ID, function(x) {
  strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])]
})
# add fake stations for climate indices
index_stations = data.frame(ID = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                   "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                   "SSTsNino4Mean.1854"),
                            longname = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                         "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                         "SSTsNino4Mean.1854"),
                            name = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                     "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                     "SSTsNino4Mean.1854"),
                            lat = -85, lon = seq(-150,150, length.out = 9),
                            ID.1 = c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                     "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                     "SSTsNino4Mean.1854"),
                            matchID = c("VAMOunsmoothedLong", "VIPOunfilteredV5", "VPDOnew", "VSAMRecon",
                                        "VSOInew", "VSSTsNino1.2Mean.1854", "VSSTsNino3.4Mean.1854", "VSSTsNino3Mean.1854",
                                        "VSSTsNino4Mean.1854"),
                            start = "1900", end = "2019", var = "ind")
all_stations = rbind(all_stations, index_stations)

# check that all names are alphabetically ordered
order(colnames(var_imp_total))
order(all_stations$matchID)

# add importance
station_importance = all_stations
# station_importance[,c("importance_King_Haakon", "importance_Ross",
#                       "importance_East_Antarctica", "importance_Weddell",
#                       "importance_Bellingshausen_Amundsen")] = t(var_imp_total)
station_importance[,paste0("M",1:nrow(var_imp_total))] = t(var_imp_total)

# aggregate by station
station_importance_agg = station_importance %>%
  group_by(name,  lat, lon) %>%
  summarise(across(M1:M24, sum))

# plot station importance
ant <- ggplot2::map_data(map = "world", region = ".",
                         # orientation = c(-90, 0, 0),
                         wrap = c(-180, 180, -90))
sectors = c("King_Hakon", "Ross",
            "East_Antarctica", "Weddell",
            "Bellingshausen_Amundsen_Sea")
pretty_sectors = c("King Haakon VII", "Ross Sea",
                   "East Antarctica", "Weddell Sea",
                   "Bellingshausen Amundsen Sea")
# first empty out plots
file.remove(list.files("plots/station_importance_fixLOOCV/sector_season", full.names = TRUE))
for (i in 1:nrow(var_imp_ss)) {
  this.sector = var_imp_ss$sector[i]
  this.pretty_sector = var_imp_ss$pretty_sector[i]
  this.season = var_imp_ss$season[i]
  this.pretty_season = c("DJF", "MAM", "JJA", "SON")[this.season]
  # make sure we get the right importance per sector
  filter_dat = station_importance_agg[,c(1:3, 3+i)]
  names(filter_dat)[4] = "importance"
  # subset to most important 25% of stations
  filter_dat = filter_dat %>% ungroup() %>% 
    filter(importance > quantile(importance, 0.75))
  # get density without the climate indices
  station_importance_agg_noInd = station_importance_agg[!(station_importance_agg$name %in%
                                                            c("AMOunsmoothedLong", "IPOunfilteredV5", "PDOnew", "SAMRecon",
                                                              "SOInew", "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854", "SSTsNino3Mean.1854",
                                                              "SSTsNino4Mean.1854")),]
  density_dat = density(station_importance_agg_noInd$lon,
                        bw = 10,
                        from = -180, to = 180,
                        weights = unlist(station_importance_agg_noInd[,3+i]) / sum(unlist(station_importance_agg_noInd[,3+i])))
  # density_dat = density(filter_dat$lon,
  #                       bw = 10,
  #                       from = -180, to = 180,
  #                       weights = unlist(filter_dat[,4]) / sum(unlist(filter_dat[,4])))
  density_dat = data.frame(x = density_dat$x, y = density_dat$y)
  # geographic projection
  my_col = rep("steelblue", 5)
  my_col[which(this.sector == sectors)] = "firebrick"
  if (this.sector == "total") {
    my_col = rep("firebrick", 5)
  }
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group),
                 data = ant,
                 fill = "gray70", colour = "gray70") +
    # coord_map("stereographic", orientation=c(-90, 0, 0), ylim = c(-90,-50)) +
    scale_x_continuous(name = "longitude", breaks = c(-180, -90, 0, 90, 180)) +
    scale_y_continuous(name = "latitude") +
    coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 0)) +
    # sea ice based sectors
    geom_vline(xintercept = c(71, 163, -110, -67, -14), lty = 2, col = "steelblue") +
    # sea ice based sectors
    geom_label(data = data.frame(x = c(117, -150, -88.5, -40.5, 28.5)[c(5,2,1,4,3)], y = -3,
                                 label = c("East\nAntarctica",
                                           "Ross Sea",
                                           "Amund.\nBel. Sea",
                                           "Weddell\nSea",
                                           "King Haakon")[c(5,2,1,4,3)]),
               aes(x = x, y = y, label = label),
               color = my_col) +
    theme_minimal() +
    # add stations
    geom_point(aes(x = lon, y = lat, size = importance, alpha = importance),
               data = filter_dat, show.legend = FALSE) +
    ggrepel::geom_text_repel(aes(x = lon, y = lat, label = name, alpha = importance),
                             max.overlaps = 20,
                             show.legend = FALSE,
                             data = filter_dat) +
    ggtitle(paste0(this.pretty_sector, ", ", this.pretty_season, ", ", "25% most important stations")) +
    scale_alpha("Variable\nImportance") +
    scale_size("Variable\nImportance") +
    geom_line(aes(x = x, y = y * 2000 - 85), data = density_dat)
  ggsave(filename = paste0("plots/station_importance_fixLOOCV/sector_season/StationImp_",
                           var_imp_ss[i, "sector"],
                           "_ssn", var_imp_ss[i, "season"], 
                           ".png"),
         width = 8, height = 5,
         device = png())
  dev.off()
}

# ##############################################################################
# how stable are the correlations over time
yrs = 1979:2020
# read out the coefficients on original level
coefs = matrix(0, nrow = length(yrs), ncol = 1 + length(all_cov_names))
colnames(coefs) = c("intercept", all_cov_names)

this.sector = "Ross"
this.n_season = 3
corr_r = data.frame(psych::corr.test(this.all_dat_s[this.all_dat_s$n_season == this.n_season, this.sector],
                                     this.all_dat_s[this.all_dat_s$n_season == this.n_season, 9:655],
                                     use = "pairwise.complete.obs")$r)
corr_p = data.frame(psych::corr.test(this.all_dat_s[this.all_dat_s$n_season == this.n_season, this.sector],
                                     this.all_dat_s[this.all_dat_s$n_season == this.n_season, 9:655],
                                     use = "pairwise.complete.obs")$r)

for (i in 1:length(yrs)) {
  this.pred_year = yrs[i]
  # delete two years before and after the year to be predicted
  del_years = -2:2 + this.pred_year
  this.all_dat_s = all_dat_s[!all_dat_s$season_year %in% del_years,]
  
  # computed variables
  # compute correlation
  # compute separate correlation matrices per season
  this.corr_test = psych::corr.test(this.all_dat_s[this.all_dat_s$n_season == this.n_season, this.sector],
                                    this.all_dat_s[this.all_dat_s$n_season == this.n_season, 9:655],
                                    use = "pairwise.complete.obs")
  # this.corr_p = this.corr_test$p
  corr_r[i,] = this.corr_test$r
  corr_p[i,] = this.corr_test$p
}
# get the correlation of 
# apply(corr_r, 2, sd) %>% hist(main = "sd(rho)")
apply(corr_p, 2, sd) %>% hist(main = "sd(p value)")

# Kerguelen in Ross Sea MAM
pdf("plots/Correlation_stability/Ross_MAM_kerguelen.pdf", height = 3, width = 4)
corr_p$V619980_SLP %>% hist(main = paste0("Kerguelen SLP, sd = ", round(sd(corr_p$V619980_SLP), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
corr_p$V619980_SLP_l1 %>% hist(main = paste0("Kerguelen SLP l1, sd = ", round(sd(corr_p$V619980_SLP_l1), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
corr_p$V619980_SLP_l2 %>% hist(main = paste0("Kerguelen SLP l2, sd = ", round(sd(corr_p$V619980_SLP_l2), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
corr_p$V619980_SLP_l3 %>% hist(main = paste0("Kerguelen SLP l3, sd = ", round(sd(corr_p$V619980_SLP_l3), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
corr_p$V619980_tmp %>% hist(main = paste0("Kerguelen TMP, sd = ", round(sd(corr_p$V619980_tmp), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
corr_p$V619980_tmp_l1 %>% hist(main = paste0("Kerguelen TMP l1, sd = ", round(sd(corr_p$V619980_tmp_l1), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
corr_p$V619980_tmp_l2 %>% hist(main = paste0("Kerguelen TMP l2, sd = ", round(sd(corr_p$V619980_tmp_l2), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
corr_p$V619980_tmp_l3 %>% hist(main = paste0("Kerguelen TMP l3, sd = ", round(sd(corr_p$V619980_tmp_l3), 2)), xlim = c(0, 1), breaks = seq(0, 1, length.out = 21))
dev.off()

corr_p$V619980_SLP %>% sd()
corr_p$V619980_SLP_l1 %>% sd()
corr_p$V619980_SLP_l2 %>% sd()
corr_p$V619980_SLP_l3 %>% sd()
corr_p$V619980_tmp %>% sd()
corr_p$V619980_tmp_l1 %>% sd()
corr_p$V619980_tmp_l2 %>% sd()
corr_p$V619980_tmp_l3 %>% sd()


for (this.sector in c("total", "King_Hakon", "Ross", "East_Antarctica",
                     "Weddell", "Bellingshausen_Amundsen_Sea")) {
  this.pretty_sector = c("total" = "Total","King_Hakon" = "King Haakon VII","Ross" = "Ross Sea", "East_Antarctica" = "East Antarctica",
                         "Weddell" = "Weddell Sea", "Bellingshausen_Amundsen_Sea" = "Bellingshausen Amundsen Sea")[this.sector]
  for (this.n_season in 1:4) {
    this.pretty_season = c("DJF", "MAM", "JJA", "SON")[this.n_season]
    for (i in 1:length(yrs)) {
      this.pred_year = yrs[i]
      # delete two years before and after the year to be predicted
      del_years = -2:2 + this.pred_year
      this.all_dat_s = all_dat_s[!all_dat_s$season_year %in% del_years,]
      
      # computed variables
      # compute correlation
      # compute separate correlation matrices per season
      this.corr_test = psych::corr.test(this.all_dat_s[this.all_dat_s$n_season == this.n_season, this.sector],
                                        this.all_dat_s[this.all_dat_s$n_season == this.n_season, 9:655],
                                        use = "pairwise.complete.obs")
      # this.corr_p = this.corr_test$p
      corr_r[i,] = this.corr_test$r
      corr_p[i,] = this.corr_test$p
    }
    png(paste0("plots/Correlation_stability/", this.sector, "_ssn", this.n_season, "sdCorPvalue.png"), 
        height = 300, width = 360)
    apply(corr_p, 2, sd) %>% hist(main = paste0(this.pretty_sector, ", ", this.pretty_season), xlim = c(0, 0.3))
    dev.off()
    png(paste0("plots/Selection_stability/", this.sector, "_ssn", this.n_season, "selProb.png"), 
        height = 300, width = 360)
    apply(corr_p < 0.05, 2, mean) %>% hist(main = paste0(this.pretty_sector, ", ", this.pretty_season), ylim = c(0, 100), xlim = c(0, 1))
    dev.off()
  }
}

# x = data.frame(prob = apply(corr_p < 0.05, 2, mean))
# ggplot(main = paste0(this.pretty_sector, ", ", this.pretty_season), xlim = c(0, 1)) +
#   geom_histogram(aes(prob), data = x, binwidth = 0.05) +
#   theme_minimal() +
#   coord_cartesian(ylim = c(0, 50)) 
# +  scale_y_log10()
