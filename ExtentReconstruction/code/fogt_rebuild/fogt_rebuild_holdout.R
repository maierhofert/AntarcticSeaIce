# this file reconstructs the analysis done by Fogt et al 2022
library("dplyr")

# read all_data_s created in fogt_rebuild.R
all_data_s = readRDS("all_data_seasonal_lagged.RDS")
all_data_s = all_data_s[all_data_s$season_year <= 2020]

# this is the holdout data
all_data_s_holdout = all_data_s[all_data_s$season_year %in% 1979:1984,]
# this is the remaining training data
all_data_s = all_data_s[all_data_s$season_year %in% 1985:2022,]
# ##############################################################################
# compute separate correlation matrices per season 
corr_test = list()
corr_p = list()
for (i in 1:4) {
  corr_test[[i]] = psych::corr.test(all_data_s[all_data_s$n_season == i, 3:8],
                                    all_data_s[all_data_s$n_season == i, 9:(ncol(all_data_s) - 1)],
                                    use = "pairwise.complete.obs")
  corr_p[[i]] = corr_test[[i]]$p
}

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

res_data$sector = as.character(res_data$sector)
res_data$pretty_sector = factor(pretty_sector, levels = pretty_sector)
# add season
res_data$season = "DJF"
res_data$season[res_data$n_season == 2] = "MAM"
res_data$season[res_data$n_season == 3] = "JJA"
res_data$season[res_data$n_season == 4] = "SON"
res_data$season = factor(res_data$season, 
                         levels = c("DJF", "MAM", "JJA", "SON"))
res_data$formula = NA
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
# it should be
sector = colnames(all_data_s)[3:8]
season = c("DJF", "MAM", "JJA", "SON")
n_season = 1:4
cor_threshold = c(0.1, 0.05, 0.025, 0.01)
n_indices = 0:9
n_pca = seq(2, 35, by = 2)

# 28.8k models
5*4*4*10*36

# fit all the models and compute validation performance
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
  this.pca = prcomp(all_dat_s[all_dat_s$season == this.season, this.covariates])
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
  
  yrs = 1985:2020
  preds = rep(NA, length(yrs))
  for (i in 1:length(yrs)) {
    this.pred_year = yrs[i]
    # delete two years before and after the year to be predicted
    del_years = -2:2 + this.pred_year
    this.mod_data_y = this.mod_data[!this.mod_data$season_year %in% del_years,]
    
    this.mod = lm(this.formula, 
                  data = this.mod_data_y)
    # mod_list[[i]] = this.mod
    # summary(this.mod)
    preds[i] = predict(this.mod, newdata = this.mod_data[this.mod_data$season_year == this.pred_year,])
  }
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

# save version with duplicates and n_pca
res_data_old = res_data

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
# saveRDS(res_data_red, "../output/fogt_rebuild_results_holdout.RDS")
saveRDS(res_data_red, "../output/fogt_rebuild_results_holdout_lindetrend.RDS")
# res_data_red = readRDS("../output/fogt_rebuild_results_holdout.RDS")


res_data_red$holdout_cor = NA
res_data_red$holdout_ce = NA
res_data_red$holdout_rmse = NA
res_data_red$holdout_sse = NA

# compute holdout performance for all unique models
for (i in 1:nrow(res_data_red)) {
  this.covariates = strsplit(res_data_red$this.covariates[i], ", ")[[1]]
  this.season = res_data_red$season[i]
  this.sector = res_data_red$sector[i]
  if (length(this.covariates) >= 2) {
    # do a PC transform
    this.pca = prcomp(all_dat_s[all_dat_s$season == this.season, this.covariates])
    
    # only retain first 35 principal components
    this.mod_data = data.frame(all_dat_s[all_dat_s$season == this.season, 
                                              1:8],
                               this.pca$x)
    this.newdata = predict(this.pca, 
                           all_data_s_holdout[all_data_s_holdout$season == this.season,]) %>%
      data.frame()
  } else if (length(this.covariates) == 1) {
    this.mod_data = all_dat_s[all_dat_s$season == this.season,]
    this.newdata = all_data_s_holdout[all_data_s_holdout$season == this.season,]
  } else {
    stop("This is not implemented")
  }
  
  this.mod = lm(res_data_red$formula[i], data = this.mod_data)
  
  preds = predict(this.mod, newdata = this.newdata)
  truth = all_data_s_holdout[all_data_s_holdout$season == this.season,this.sector]
  res_data_red[i,"holdout_cor"] = cor(preds, truth)
  res_data_red[i,"holdout_ce"] = my_CE(truth, preds)
  res_data_red[i,"holdout_rmse"] = sqrt(mean((preds - truth)^2))
  res_data_red[i,"holdout_sse"] = sum((preds - truth)^2)
}
summary(res_data_red)

# for some models this just goes horribly wrong
test = res_data_red[res_data_red$holdout_sse > 100,]

plot(res_data_red$val_ce, res_data_red$holdout_ce, ylim = c(-1, 1), xlim = c(0.5, 1))
plot(res_data_red$val_cor, res_data_red$holdout_cor, xlim = c(0.5, 1))

ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  # facet_grid(season~pretty_sector) +
  geom_point(aes(x = val_cor, y = holdout_cor)) +
  scale_x_continuous("Validation Correlation") + 
  scale_y_continuous("Holdout Correlation") +
  theme_bw()
ggsave("plots/fogt_rebuild_holdoutCor.png", height = 4, width = 7)
ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_grid(season~pretty_sector) +
  geom_point(aes(x = val_cor, y = holdout_cor)) +
  scale_x_continuous("Validation Correlation") + 
  scale_y_continuous("Holdout Correlation") +
  theme_bw()
ggsave("plots/fogt_rebuild_holdoutCor_grid.png", height = 7, width = 12)


ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  # facet_grid(season~pretty_sector) +
  geom_point(aes(x = val_ce, y = holdout_ce)) +
  scale_x_continuous("Validation CE") + 
  scale_y_continuous("Holdout CE", limits = c(-1, 1)) +
  theme_bw()
ggsave("plots/fogt_rebuild_holdoutCE.png", height = 4, width = 7)
ggplot(res_data_red[res_data_red$sector != "total",]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_grid(season~pretty_sector) +
  geom_point(aes(x = val_ce, y = holdout_ce)) +
  scale_x_continuous("Validation CE") + 
  scale_y_continuous("Holdout CE", limits = c(-1, 1)) +
  theme_bw()
ggsave("plots/fogt_rebuild_holdoutCE_grid.png", height = 7, width = 12)

# how many models do we have per sector and season
res_data_red %>%
  group_by(season, sector) %>%
  summarise(n())
tab = table(res_data_red[,c("sector", "season")])
tab
colSums(tab[-5,])

# print out more summary statistics
res_data_red %>%
  group_by(season, sector) %>%
  summarise(holdout_ce = mean(holdout_ce))

res_data_red %>%
  group_by(sector) %>%
  summarise(holdout = mean(holdout_ce))

res_data_red %>%
  group_by(sector) %>%
  summarise(holdout_cor = mean(holdout_cor))

res_data_red %>%
  filter(sector != "Total") %>%
  group_by(sector) %>%
  summarise(holdout_cor = median(holdout_cor)) %>%
  summarise(holdout_cor = median(holdout_cor))
# 0.539

res_data_red %>%
  filter(sector != "Total") %>%
  group_by(sector) %>%
  summarise(holdout_ce = median(holdout_ce)) %>%
  summarise(holdout_ce = median(holdout_ce))
# -0.733
