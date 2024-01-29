# this file reconstructs the analysis done by Fogt et al 2022
library("dplyr")
library("ggplot2")

# read monthly data
all_data = readRDS("all_data.RDS")
names(all_data)[-(1:9)] = paste0("V", names(all_data)[-(1:9)])
names(all_data)
#  long record 
long_record = c("619010_SLP", "685880_SLP", "688160_SLP", "688420_SLP",
                "837430_SLP", "855740_SLP", "855850_SLP", "857660_SLP",
                "859340_SLP", "870470_SLP", "872220_SLP", "873440_SLP",
                "875850_SLP", "877500_SLP", "878490_SLP", "889680_SLP",
                "919380_SLP", "931190_SLP", "934340_SLP", "936150_SLP",
                "937800_SLP", "938940_SLP", "939870_SLP", "943260_SLP",
                "945780_SLP", "946100_SLP", "946720_SLP", "947680_SLP",
                "948680_SLP", "949750_SLP",
                "619010_tmp", "685880_tmp", "688160_tmp", "688420_tmp",
                "837430_tmp", "855740_tmp", "855850_tmp", "857660_tmp",
                "859340_tmp", "870470_tmp", "872220_tmp", "873440_tmp",
                "875850_tmp", "877500_tmp", "878490_tmp", "889680_tmp",
                "919380_tmp", "931190_tmp", "934340_tmp", "936150_tmp",
                "937800_tmp", "938940_tmp", "939870_tmp", "943260_tmp",
                "945780_tmp", "946100_tmp", "946720_tmp", "947680_tmp",
                "948680_tmp", "949750_tmp",
                "IPOunfilteredV5", "AMOunsmoothedLong", "PDOnew", "SOInew",
                "SSTsNino1.2Mean.1854", "SSTsNino3.4Mean.1854",
                "SSTsNino3Mean.1854", "SSTsNino4Mean.1854", "SAMRecon")
long_record_V = paste0("V", long_record)
# subset to long record variables
all_data = all_data[, c(1:9, which(names(all_data) %in% long_record_V))]
# which(!(long_record_V %in% names(all_data)))
# long_record_V[47]

# what happens if we detrend all covariates and sea ice by sector
all_data = all_data[which(!is.na(all_data$total)),]
slope = rep(NA, ncol(all_data))
for (i in c(10:ncol(all_data))) { # c(3:8, 10:173) # do not include sea ice in the standardization
  this.lm = lm(all_data[,i] ~ 1 + I(1:nrow(all_data)))
  slope[i] = coef(this.lm)[2]
  all_data[!is.na(all_data[,i]),i] = residuals(this.lm)
}
# there is not much going on in terms of standardization
hist(slope)


# add season
all_data$n_season = 1
all_data$n_season[all_data$Month %in% 3:5] = 2
all_data$n_season[all_data$Month %in% 6:8] = 3
all_data$n_season[all_data$Month %in% 9:11] = 4

# make sure December gets combined with next years January and February
all_data$season_year = lead(all_data$Year, 1)
all_data$season_year[nrow(all_data)] = 2022

# ####################################################################
# add lagged versions of the variables
cov = all_data[,-c(1:9)]
lags = c(1, 2, 3)
col_names = names(cov)
jmax = ncol(cov)
# add buffer for lagging
cov = rbind(matrix(rnorm(12 * jmax),
                   ncol = jmax, nrow = 12,
                   dimnames = list(NULL, colnames(cov))),
            cov)
for (j in 1:jmax) { # iterate over all covariates
  for (i in seq_along(lags)) { # for every covariate, iterate over all lags
    cov[,paste0(col_names[j], "_l", lags[i])] = lag(cov[,col_names[j]], n = lags[i])
  }
}
cov = cov[-(1:12),]
cov$date = all_data$Date
all_data = right_join(all_data[,1:9],
                      cov,
                      by = c(Date = "date"))
# remove missing half of 2022
all_data = all_data[!(all_data$Year == 2022 & all_data$Month >= 6),]

# check what we have done with the lags for one variable
# all_data[1452:1464, c("Year", "Month",  "season_year", "n_season", "V619010_SLP", "V619010_SLP_l1", "V619010_SLP_l2", "V619010_SLP_l3")] %>% 
#   round(2)
all_data[495:507, c("Year", "Month",  "season_year", "n_season", "V619010_SLP", "V619010_SLP_l1", "V619010_SLP_l2", "V619010_SLP_l3")] %>% 
  round(2)
# ##########################################################################
# compute seasonal mean
all_data_s = all_data %>%
  group_by(season_year, n_season) %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), mean)) %>%
  select(-Year, -Month, 
         -n_season_l1, -n_season_l2, -n_season_l3,
         -season_year_l1, -season_year_l2, -season_year_l3)
# get labeled seasons
all_data_s$season = "DJF"
all_data_s$season[all_data_s$n_season == 2] = "MAM"
all_data_s$season[all_data_s$n_season == 3] = "JJA"
all_data_s$season[all_data_s$n_season == 4] = "SON"

# save seasonal data
all_data_s = ungroup(all_data_s)
all_data_s = as.data.frame(all_data_s)

# compare 2020 in seasonal data to 2020 in monthly data
# all_data[1452:1464, c("Year", "Month",  "season_year", "n_season", "V619010_SLP", "V619010_SLP_l1", "V619010_SLP_l2", "V619010_SLP_l3")] %>% 
#   round(2)
# all_data_s[485:488, c("season_year", "n_season", "V619010_SLP", "V619010_SLP_l1", "V619010_SLP_l2", "V619010_SLP_l3")] %>% 
#   round(2)
all_data[495:507, c("Year", "Month",  "season_year", "n_season", "V619010_SLP", "V619010_SLP_l1", "V619010_SLP_l2", "V619010_SLP_l3")] %>% 
  round(2)
all_data_s[166:169, c("season_year", "n_season", "V619010_SLP", "V619010_SLP_l1", "V619010_SLP_l2", "V619010_SLP_l3")] %>% 
  round(2)

# saveRDS(all_data_s, "../data/complete_dataset/all_data_seasonal_lagged.RDS")
# saveRDS(all_data_s, "all_data_seasonal_lagged.RDS")
saveRDS(all_data_s, "all_data_seasonal_lagged_longRecord.RDS")
# subset data to 1979 to 2020
all_dat_s = all_data_s[between(all_data_s$season_year, 1979, 2020),]

# # permute sea ice
# all_dat_s[,3:8] = all_dat_s[sample(1:nrow(all_dat_s), nrow(all_dat_s)),3:8]

# ##############################################################################
# compute separate correlation matrices and p values per season 
corr_test = list()
corr_p = list()
for (i in 1:4) {
  corr_test[[i]] = psych::corr.test(all_data_s[all_data_s$n_season == i, 3:8],
                                    all_data_s[all_data_s$n_season == i, 9:(ncol(all_data_s) - 1)],
                                    use = "pairwise.complete.obs")
  corr_p[[i]] = corr_test[[i]]$p
}
# look into the actual covariance matrices
# look for DJF (n_season = 1), Bellingshausen Amundsen Sea (6th sector including total)
# station_data = read.csv("../../data/fogt_predictor_data/station_data.csv")

# with station DURBAN/ LOUIS BOT 685880_SLP 
corr_test[[1]]$r[6,c("V685880_SLP", "V685880_SLP_l1", "V685880_SLP_l2", "V685880_SLP_l3")]

# Dunedin
corr_test[[1]]$r[6,c("V938940_SLP", "V938940_SLP_l1", "V938940_SLP_l2", "V938940_SLP_l3")]


# helper function, coefficient of efficiency
my_CE = function(truth, pred) {
  1 - sum((truth - pred)^2) / sum((truth - mean(truth))^2)
}

# ##############################################################################
# set up a nice cross validation
# sectors
sector = colnames(all_data_s)[3:8]
pretty_sector = c("Total", "King Haakon VII", "Ross Sea", "East Antarctica", 
                  "Weddell Sea", "Bellingshausen Amundsen Sea")
# season
season = c("DJF", "MAM", "JJA", "SON")
n_season = 1:4
# correlation threshold
cor_threshold = c(0.1, 0.05, 0.025, 0.01)
# number of climate indices used
n_indices = 0:9
# only fit every other number of pc retained to speed up computation
n_pca = seq(2, 35, by = 2)


# set up cross validation results data
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

# make results data pretty
res_data$sector = as.character(res_data$sector)
res_data$pretty_sector = factor(pretty_sector, levels = pretty_sector)
# add season
res_data$season = "DJF"
res_data$season[res_data$n_season == 2] = "MAM"
res_data$season[res_data$n_season == 3] = "JJA"
res_data$season[res_data$n_season == 4] = "SON"
res_data$season = factor(res_data$season, 
                         levels = c("DJF", "MAM", "JJA", "SON"))
# formula for linear regression model
res_data$formula = NA

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

# add covariate prediction strength
all_cov_names = names(all_dat_s)[9:(ncol(all_dat_s) - 1)]
res_data$intercept = 0
res_data[,all_cov_names] = 0

# it should be 28.8k models, but we only do every other number of PC retained
5*4*4*10*36

# actual cross validation, iterate through results data and fit each model
for (j in 1:nrow(res_data)) {
  # print progress
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
  
  
  # computed variables
  # which covariates meet the correlation significance threshold in this season
  this.covariates = colnames(corr_p[[this.n_season]])[corr_p[[this.n_season]][this.sector,] < this.cor_threshold]
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
  
  # do a PC transform
  this.pca = prcomp(all_dat_s[all_dat_s$season == this.season, this.covariates], center = FALSE)
  # which PCs are the most strongly correlated?
  this.pc_cor = cor(all_dat_s[all_dat_s$season == this.season, this.sector],
                    this.pca$x[,1:min(length(this.covariates), 35)])
  # order the PC by correlation magnitude
  this.strongest_pc = colnames(this.pc_cor)[order(abs(this.pc_cor), decreasing = TRUE)]
  # get the strongest corelations
  this.strongest_pc = this.strongest_pc[1:min(this.n_pca, length(this.covariates))]
  
  # only retain first 35 principal components
  this.mod_data = data.frame(all_dat_s[all_dat_s$season == this.season, 
                                       1:8],
                             this.pca$x[,this.strongest_pc])
  
  # create formula
  this.formula = paste0(this.sector,
                        " ~ ",
                        paste(this.strongest_pc,
                              collapse = " + ", sep = "")) %>%
    formula()
  
  yrs = 1979:2020
  preds = rep(NA, length(yrs))
  my_pred = rep(NA, length(yrs))
  # read out the coefficients on original level
  coefs = matrix(NA, nrow = length(yrs), ncol = 1 + length(this.covariates))
  # all_data_s[all_data_s$n_season == i, 9:656]
  for (i in 1:length(yrs)) {
    this.pred_year = yrs[i]
    # delete two years before and after the year to be predicted
    del_years = -2:2 + this.pred_year
    this.mod_data_y = this.mod_data[!this.mod_data$season_year %in% del_years,]
    
    this.mod = lm(this.formula, 
                  data = this.mod_data_y)
    # # confirm T = X W
    # this.pca$x - as.matrix(all_dat_s[all_dat_s$season == this.season, this.covariates]) %*% this.pca$rotation
    # this.pca$x - predict(this.pca, newdata = all_dat_s[all_dat_s$season == this.season, this.covariates])
    # this.pca$x
    
    # mod_list[[i]] = this.mod
    # summary(this.mod)
    preds[i] = predict(this.mod, newdata = this.mod_data[this.mod_data$season_year == this.pred_year,])
    
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
    coefs[i,] = c(this.coef[1], new_coef)
  }
  # compute average strength of covariates and add it to the results
  res_data[j, c("intercept", this.covariates)] = colMeans(abs(coefs))

  # # plot predicted vs fitted
  # plot(1979:2020, this.mod_data[,this.sector], type = "b", col = "green")
  # abline(h = 0, lty = 2)
  # points(1979:2020, preds, type = "b")
  
  res_data[j,"this.covariates"] = paste(this.covariates, collapse = ", ")
  res_data[j,"formula"] = deparse(this.formula, width.cutoff = 500)
  res_data[j,"val_cor"] = cor(preds, this.mod_data[,this.sector])
  res_data[j,"val_ce"] = my_CE(this.mod_data[,this.sector], preds)
  res_data[j,"val_rmse"] = sqrt(mean((preds - this.mod_data[,this.sector])^2))
  res_data[j,"rmse0"] = sqrt(mean((this.mod_data[,this.sector])^2))
}

# save version with duplicates and n_pca in object res_data_old
res_data_old = res_data
# res_data = res_data_old

# # order by n_pca performance
# res_data2 = res_data %>%
#   arrange(sector, n_season, cor_threshold, n_indices, desc(val_ce))

# find the optimal number of principal components to be retained
res_data3 = res_data %>%
  group_by(sector, n_season, cor_threshold, n_indices) %>%
  slice(which.max(val_ce))
table(res_data3$n_pca)
duplicated(res_data3[, c("formula", "n_season")]) %>% sum()

# remove duplicated models
res_data4 = res_data3[!duplicated(res_data3[, c("formula", "n_season")]),]
table(res_data4$n_pca)
hist(res_data4$n_pca)

# reduced data to only unique models
res_data_red = res_data4
hist(res_data_red$val_ce[res_data_red$val_ce > 0])
hist(res_data_red$val_cor)

# save without duplicates
saveRDS(res_data_red, "fogt_rebuild_results_lindetrend_permute4.RDS")
# res_data_red = readRDS("fogt_rebuild_results_lindetrend.RDS")

library("ggplot2")

# plot validation correlation
ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~season, nrow = 2) +
  geom_boxplot(aes(x = pretty_sector, y = val_cor)) +
  scale_x_discrete("Sector") + 
  scale_y_continuous("Validation Correlation", limits = c(-0.5, 1)) +
  theme_bw()
ggsave("plots/fogt_rebuild_valCor.png", height = 7, width = 12)

# plot validation CE
ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~season, nrow = 2) +
  geom_boxplot(aes(x = pretty_sector, y = val_ce)) +
  scale_x_discrete("Sector") +
  scale_y_continuous("Validation CE", limits = c(-0.5, 1)) +
  theme_bw()
ggsave("plots/fogt_rebuild_valCe.png", height = 7, width = 12)

# how many models do we have per sector and season
res_data_red %>%
  group_by(season, sector) %>%
  summarise(n())
tab = table(res_data_red[,c("sector", "season")])
tab
# number of models per season across all sectors
colSums(tab[-5,])

# average validation CE by season and sector
res_data_red %>%
  group_by(season, sector) %>%
  summarise(val_ce = mean(val_ce))

# average validation CE by sector
res_data_red %>%
  group_by(sector) %>%
  summarise(val_ce = mean(val_ce))

# average validation correlation by sector
res_data_red %>%
  group_by(sector) %>%
  summarise(val_cor = mean(val_cor))

# overall average validation correlation
res_data_red %>%
  filter(sector != "Total") %>%
  group_by(sector) %>%
  summarise(val_cor = mean(val_cor)) %>%
  summarise(val_cor = mean(val_cor))
# 0.843

# overall average validation CE
res_data_red %>%
  filter(sector != "Total") %>%
  group_by(sector) %>%
  summarise(val_ce = mean(val_ce)) %>%
  summarise(val_ce = mean(val_ce))
# 0.700


# ##############################################################################
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
file.remove(list.files("plots/station_importance", full.names = TRUE))
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
  ggsave(filename = paste0("plots/station_importance/StationImp_",
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
file.remove(list.files("plots/station_importance/sector_season", full.names = TRUE))
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
  ggsave(filename = paste0("plots/station_importance/sector_season/StationImp_",
                           var_imp_ss[i, "sector"],
                           "_ssn", var_imp_ss[i, "season"], 
                           ".png"),
         width = 8, height = 5,
         device = png())
  dev.off()
}



