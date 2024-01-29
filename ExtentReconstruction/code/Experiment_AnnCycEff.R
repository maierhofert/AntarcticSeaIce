# this file fits frequentist models to the data
# combine covariates and data
source("target_sea_ice.R")
source("cov_aggregation.R")

all_data = full_join(cbind(data.frame(tdate, Doy, Date = sie_date),
                           sie_detrended),
                     cbind(cov, 
                           eigen_x,
                           by = "Date"))
all_data = all_data[order(all_data$Date),]

all_data$Date[!is.na(all_data$King_Hakon)] %>% head()

# subset data to after 1989
all_dat = all_data[all_data$Date >= as.Date("1979-01-01"),]
# all_dat$Doy = round(all_dat$Doy)

# # plot principal components to see how smooth they are
# plot(all_dat$Date[test_id], all_dat$PC1[test_id], type = "l")
# plot(all_dat$Date[test_id], all_dat$PC2[test_id], type = "l")
# plot(all_dat$Date[test_id], all_dat$PC5[test_id], type = "l")
# plot(all_dat$Date[test_id], all_dat$PC10[test_id], type = "l")
# plot(all_dat$Date[test_id], all_dat$PC20[test_id], type = "l")
# plot(all_dat$Date[test_id], all_dat$PC50[test_id], type = "l")

# set up a nice cross validation
# sectors
sector = colnames(sie_detrended)
pretty_sector = c("King Hakon VII", "Ross Sea", "East Antarctica", "Weddell Sea", 
                  "Bellingshausen\nAmundsen Sea")
# blocks for cross validation
year_block = list(1979:1984, 1985:1988,
                  1989:1994, 1995:1998,
                  1999:2004, 2005:2009,
                  2010:2014, 2015:2018)
year_block_ind = 1:length(year_block)
# number of principal components to use
n_pc = seq(1, 31, by = 5)

# set up results data
res_data = expand.grid(sector = sector, year_block_ind = year_block_ind, n_pc = n_pc)
res_data$pretty_sector = pretty_sector
# compute rmse
res_data$rmse = NA
res_data$rmse0 = NA
res_data$cor = NA
res_data$rmse_train = NA
rmse = function(x1, x2) {
  sqrt(mean((x1 - x2)^2, na.rm = TRUE))
}
mod_list = list()
for (i in 1:nrow(res_data)) {
  # split into training and test data
  # get firs and last year of test data
  first_year = year_block[[res_data[i,"year_block_ind"]]] %>% min()
  last_year = year_block[[res_data[i,"year_block_ind"]]] %>% max()
  # get indices for test and training data
  test_id = between(all_dat$Date, 
                    as.Date(paste0(first_year, "-01-01")), 
                    as.Date(paste0(last_year,"-01-01")))
  train_id = !test_id
  # create formula
  this.sector = res_data[i,"sector"] %>% as.character()
  this.pretty_sector = res_data[i,"pretty_sector"] %>% as.character()
  # this.formula = paste0(this.sector, 
  #                       " ~ ", 
  #                       paste(
  #                         paste0("s(Doy, by = PC", 
  #                                1:res_data[i,"n_pc"],
  #                                ", bs = 'cc')"), collapse = " + ")) %>%
  #   formula()
  # fit model with linear effect per day of year
  # lm_ab_annualcycle = mgcv::gam(this.formula,
  #                               data = all_dat[train_id,])
  
  # # chuck it into random forest instead
  this.formula = paste0(this.sector, 
                        " ~ Doy + ", 
                        paste(
                          paste0("PC", 
                                 1:res_data[i,"n_pc"]), collapse = " + ")) %>%
    formula()
  lm_ab_annualcycle = randomForest::randomForest(this.formula,
                                                 sampsize = 5000,
                                                 nodesize = 15,
                                data = all_dat[train_id, ], na.action = na.omit)
  mod_list[[i]] = lm_ab_annualcycle
  # plot(lm_ab_annualcycle, scheme = 2, se = F)
  # anova(lm_ab_annualcycle)
  
  preds = predict(lm_ab_annualcycle,
                  all_dat[test_id,])
  res_data$cor[i] = cor(preds,
                        all_dat[test_id, this.sector], use = "pairwise")
  res_data$rmse0[i] = rmse(0, all_dat[test_id, this.sector])
  res_data$rmse[i] = rmse(preds, all_dat[test_id, this.sector])
  
  res_data$rmse_train[i] = rmse(predict(lm_ab_annualcycle,
                                        all_dat[train_id,]), 
                                all_dat[train_id, this.sector])
  
  # plot(all_dat$Date[test_id],
  #      all_dat[test_id, this.sector], type = "l",
  #      ylab = "sea ice anomaly", xlab = "Date",
  #      main = this.pretty_sector)
  # lines(all_dat$Date[test_id], preds, col = "blue")
  # abline(h = 0, lty = 2)
}
# saveRDS(mod_list, "../data/models/mod_list_AnnCycLiMo.RDS")
# saveRDS(res_data, "../data/models/res_data_AnnCycLiMo.RDS")
saveRDS(mod_list, "../data/models/mod_list_RF.RDS")
saveRDS(res_data, "../data/models/res_data_RF.RDS")

# plot results
ggplot(res_data) +
  facet_wrap(vars(pretty_sector)) +
  geom_boxplot(aes(x = factor(n_pc), y = rmse_train, col = pretty_sector)) +
  scale_color_discrete("Sector") +
  scale_x_discrete("Number of PCs used") +
  scale_y_continuous("training RMSE") +
  theme_minimal()
# ggsave("../plots/performance/AnnCycLiMo_trainrmse.png", device = "png", width = 7, height = 4)
ggsave("../plots/performance/RF_trainrmse.png", device = "png", width = 7, height = 4)

ggplot(res_data) +
  facet_wrap(vars(pretty_sector)) +
  geom_boxplot(aes(x = factor(n_pc), y = rmse, col = pretty_sector)) +
  scale_color_discrete("Sector") +
  scale_x_discrete("Number of PCs used") +
  scale_y_continuous("test RMSE") +
  theme_minimal()
#ggsave("../plots/performance/AnnCycLiMo_testrmse.png", device = "png", width = 7, height = 4)
ggsave("../plots/performance/RF_testrmse.png", device = "png", width = 7, height = 4)


ggplot(res_data) +
  facet_wrap(vars(pretty_sector)) +
  geom_boxplot(aes(x = factor(n_pc), y = rmse0 - rmse, col = pretty_sector)) +
  scale_color_discrete("Sector") +
  scale_x_discrete("Number of PCs used") +
  scale_y_continuous("predict 0 RMSE - test RMSE") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal()
# ggsave("../plots/performance/AnnCycLiMo_rmsediff.png", device = "png", width = 7, height = 4)
ggsave("../plots/performance/RF_rmsediff.png", device = "png", width = 7, height = 4)

# check out variable importance
varImpPlot(mod_list[[276]], main = res_data$pretty_sector[276])
varImpPlot(mod_list[[277]], main = res_data$pretty_sector[277])
varImpPlot(mod_list[[278]], main = res_data$pretty_sector[278])
varImpPlot(mod_list[[279]], main = res_data$pretty_sector[279])
varImpPlot(mod_list[[280]], main = res_data$pretty_sector[280])



# ##############################################################################
# split into training and test data
test_id = between(all_dat$Date, as.Date("1989-01-01"), as.Date("1994-01-01"))
train_id = !test_id

lm_ab_annualcycle2 = mgcv::gam(King_Hakon ~ 
                                 # can I still include an intercept?
                                 # I know I detrended the data
                                 # s(Doy, bs = "cc") +
                                 s(Doy, by = PC1, bs = "cc") +
                                 s(Doy, by = PC2, bs = "cc") +
                                 s(Doy, by = PC3, bs = "cc") +
                                 s(Doy, by = PC4, bs = "cc") +
                                 s(Doy, by = PC5, bs = "cc") +
                                 s(Doy, by = PC6, bs = "cc") +
                                 s(Doy, by = PC7, bs = "cc") +
                                 s(Doy, by = PC8, bs = "cc") +
                                 s(Doy, by = PC9, bs = "cc") +
                                 s(Doy, by = PC10, bs = "cc") +
                                 s(Doy, by = PC11, bs = "cc") +
                                 s(Doy, by = PC12, bs = "cc"), # +
                               # s(Doy, by = PC13, bs = "cc") +
                               # s(Doy, by = PC14, bs = "cc") +
                               # s(Doy, by = PC15, bs = "cc") +
                               # s(Doy, by = PC16, bs = "cc") +
                               # s(Doy, by = PC17, bs = "cc") +
                               # s(Doy, by = PC18, bs = "cc") +
                               # s(Doy, by = PC19, bs = "cc") +
                               # s(Doy, by = PC20, bs = "cc") +
                               # s(Doy, by = PC21, bs = "cc") +
                               # s(Doy, by = PC22, bs = "cc"),
                               data = all_dat[train_id,])
# plot(lm_ab_annualcycle2, scheme = 2, se = F)
anova(lm_ab_annualcycle2)

rmse = function(x1, x2) {
  sqrt(mean((x1 - x2)^2, na.rm = TRUE))
}
plot(all_dat$Date[test_id], 
     all_dat[test_id, "King_Hakon"], type = "l")
preds2 = predict(lm_ab_annualcycle2, 
                 all_dat[test_id,])
cor(preds2,
    all_dat[test_id, "King_Hakon"], use = "pairwise")
lines(all_dat$Date[test_id], preds2, col = 3)
abline(h = 0, lty = 2)
rmse(0, all_dat[test_id, "King_Hakon"])
rmse(preds2, all_dat[test_id, "King_Hakon"])

# Ross
lm_ab_annualcycle2 = mgcv::gam(Ross ~ s(Doy, by = PC1, bs = "cc") +
                                 s(Doy, by = PC2, bs = "cc") +
                                 s(Doy, by = PC3, bs = "cc") +
                                 s(Doy, by = PC4, bs = "cc") +
                                 s(Doy, by = PC5, bs = "cc") +
                                 s(Doy, by = PC6, bs = "cc") +
                                 s(Doy, by = PC7, bs = "cc") +
                                 s(Doy, by = PC8, bs = "cc") +
                                 s(Doy, by = PC9, bs = "cc") +
                                 s(Doy, by = PC10, bs = "cc") +
                                 s(Doy, by = PC11, bs = "cc") +
                                 s(Doy, by = PC12, bs = "cc"), #  +
                               # s(Doy, by = PC13, bs = "cc") +
                               # s(Doy, by = PC14, bs = "cc") +
                               # s(Doy, by = PC15, bs = "cc") +
                               # s(Doy, by = PC16, bs = "cc") +
                               # s(Doy, by = PC17, bs = "cc") +
                               # s(Doy, by = PC18, bs = "cc") +
                               # s(Doy, by = PC19, bs = "cc") +
                               # s(Doy, by = PC20, bs = "cc") +
                               # s(Doy, by = PC21, bs = "cc") +
                               # s(Doy, by = PC22, bs = "cc"), 
                               data = all_dat[train_id,])
# plot(lm_ab_annualcycle2, scheme = 2, se = F)
anova(lm_ab_annualcycle2)
preds2 = predict(lm_ab_annualcycle2, 
                 all_dat[test_id,])
cor(preds2,
    all_dat[test_id, "Ross"], use = "pairwise")
plot(all_dat$Date[test_id], 
     all_dat[test_id, "Ross"], type = "l")
lines(all_dat$Date[test_id], preds2, col = 3)
abline(h = 0, lty = 2)
rmse(0, all_dat[test_id, "Ross"])
rmse(preds, all_dat[test_id, "Ross"])
rmse(preds2, all_dat[test_id, "Ross"])

lm_ab_annualcycle2 = mgcv::gam(East_Antarctica ~ s(Doy, by = PC1, bs = "cc") +
                                 s(Doy, by = PC2, bs = "cc") +
                                 s(Doy, by = PC3, bs = "cc") +
                                 s(Doy, by = PC4, bs = "cc") +
                                 s(Doy, by = PC5, bs = "cc") +
                                 s(Doy, by = PC6, bs = "cc") +
                                 s(Doy, by = PC7, bs = "cc") +
                                 s(Doy, by = PC8, bs = "cc") +
                                 s(Doy, by = PC9, bs = "cc") +
                                 s(Doy, by = PC10, bs = "cc") +
                                 s(Doy, by = PC11, bs = "cc") +
                                 s(Doy, by = PC12, bs = "cc"), # +
                               # s(Doy, by = PC13, bs = "cc") +
                               # s(Doy, by = PC14, bs = "cc") +
                               # s(Doy, by = PC15, bs = "cc") +
                               # s(Doy, by = PC16, bs = "cc") +
                               # s(Doy, by = PC17, bs = "cc") +
                               # s(Doy, by = PC18, bs = "cc") +
                               # s(Doy, by = PC19, bs = "cc") +
                               # s(Doy, by = PC20, bs = "cc") +
                               # s(Doy, by = PC21, bs = "cc") +
                               # s(Doy, by = PC22, bs = "cc"),
                               data = all_dat[train_id,])
# plot(lm_ab_annualcycle2, scheme = 2, se = F)
anova(lm_ab_annualcycle2)
preds2 = predict(lm_ab_annualcycle2, 
                 all_dat[test_id,])
cor(preds2,
    all_dat[test_id, "East_Antarctica"], use = "pairwise")

plot(all_dat$Date[test_id], 
     all_dat[test_id, "East_Antarctica"], type = "l")
lines(all_dat$Date[test_id], preds2, col = 3)
abline(h = 0, lty = 2)
rmse(0, all_dat[test_id, "East_Antarctica"])
rmse(preds2, all_dat[test_id, "East_Antarctica"])

lm_ab_annualcycle2 = mgcv::gam(Weddell ~  s(Doy, by = PC1, bs = "cc") +
                                 s(Doy, by = PC2, bs = "cc") +
                                 s(Doy, by = PC3, bs = "cc") +
                                 s(Doy, by = PC4, bs = "cc") +
                                 s(Doy, by = PC5, bs = "cc") +
                                 s(Doy, by = PC6, bs = "cc") +
                                 s(Doy, by = PC7, bs = "cc") +
                                 s(Doy, by = PC8, bs = "cc") +
                                 s(Doy, by = PC9, bs = "cc") +
                                 s(Doy, by = PC10, bs = "cc") +
                                 s(Doy, by = PC11, bs = "cc") +
                                 s(Doy, by = PC12, bs = "cc"), # +
                               # s(Doy, by = PC13, bs = "cc") +
                               # s(Doy, by = PC14, bs = "cc") +
                               # s(Doy, by = PC15, bs = "cc") +
                               # s(Doy, by = PC16, bs = "cc") +
                               # s(Doy, by = PC17, bs = "cc") +
                               # s(Doy, by = PC18, bs = "cc") +
                               # s(Doy, by = PC19, bs = "cc") +
                               # s(Doy, by = PC20, bs = "cc") +
                               # s(Doy, by = PC21, bs = "cc") +
                               # s(Doy, by = PC22, bs = "cc"), 
                               data = all_dat[train_id,])
# plot(lm_ab_annualcycle2, scheme = 2, se = F)

preds2 = predict(lm_ab_annualcycle2, 
                 all_dat[test_id,])
cor(preds2,
    all_dat[test_id, "Weddell"], use = "pairwise")
plot(all_dat$Date[test_id], 
     all_dat[test_id, "Weddell"], type = "l")
lines(all_dat$Date[test_id], preds2, col = 3)
abline(h = 0, lty = 2)
rmse(0, all_dat[test_id, "Weddell"])
rmse(preds, all_dat[test_id, "Weddell"])
rmse(preds2, all_dat[test_id, "Weddell"])

lm_ab_annualcycle2 = mgcv::gam(Bellingshausen_Amundsen_Sea ~ s(Doy, by = PC1, bs = "cc") +
                                 s(Doy, by = PC2, bs = "cc") +
                                 s(Doy, by = PC3, bs = "cc") +
                                 s(Doy, by = PC4, bs = "cc") +
                                 s(Doy, by = PC5, bs = "cc") +
                                 s(Doy, by = PC6, bs = "cc") +
                                 s(Doy, by = PC7, bs = "cc") +
                                 s(Doy, by = PC8, bs = "cc") +
                                 s(Doy, by = PC9, bs = "cc") +
                                 s(Doy, by = PC10, bs = "cc") +
                                 s(Doy, by = PC11, bs = "cc") +
                                 s(Doy, by = PC12, bs = "cc") +
                                 s(Doy, by = PC13, bs = "cc") +
                                 s(Doy, by = PC14, bs = "cc") +
                                 s(Doy, by = PC15, bs = "cc") +
                                 s(Doy, by = PC16, bs = "cc") +
                                 s(Doy, by = PC17, bs = "cc") +
                                 s(Doy, by = PC18, bs = "cc") +
                                 s(Doy, by = PC19, bs = "cc") +
                                 s(Doy, by = PC20, bs = "cc"),
                               # s(Doy, by = PC21, bs = "cc") +
                               # s(Doy, by = PC22, bs = "cc"), 
                               data = all_dat[train_id,])
# plot(lm_ab_annualcycle2, scheme = 2, se = F)
preds2 = predict(lm_ab_annualcycle2, 
                 all_dat[test_id,])
cor(preds2,
    all_dat[test_id, "Bellingshausen_Amundsen_Sea"], use = "pairwise")
plot(all_dat$Date[test_id], 
     all_dat[test_id, "Bellingshausen_Amundsen_Sea"], type = "l")
lines(all_dat$Date[test_id], preds2, col = 3)
abline(h = 0, lty = 2)
rmse(0, all_dat[test_id, "Bellingshausen_Amundsen_Sea"])
rmse(preds, all_dat[test_id, "Bellingshausen_Amundsen_Sea"])
rmse(preds2, all_dat[test_id, "Bellingshausen_Amundsen_Sea"])













# # lm_ab_annualcycle = mgcv::gam(King_Hakon ~ s(PC1 ~ Doy, bs = "cc"), data = all_dat)
# lm_ab_annualcycle = mgcv::gam(King_Hakon ~ te(Doy, PC1, bs = c("cc", "ds")) +
#                                 te(Doy, PC2, bs = c("cc", "ds")) +
#                                 te(Doy, PC3, bs = c("cc", "ds")) +
#                                 te(Doy, PC4, bs = c("cc", "ds")) +
#                                 te(Doy, PC5, bs = c("cc", "ds")),
#                               data = all_dat[train_id,])
# # plot(lm_ab_annualcycle, scheme = 2, se = F)
# anova(lm_ab_annualcycle)
# 
# lm_ab_annualcycle2 = mgcv::gam(King_Hakon ~ 
#                                  # can I still include an intercept?
#                                  # I know I detrended the data
#                                  # s(Doy, bs = "cc") +
#                                  s(Doy, by = PC1, bs = "cc") +
#                                  s(Doy, by = PC2, bs = "cc") +
#                                  s(Doy, by = PC3, bs = "cc") +
#                                  s(Doy, by = PC4, bs = "cc") +
#                                  s(Doy, by = PC5, bs = "cc") +
#                                  s(Doy, by = PC6, bs = "cc") +
#                                  s(Doy, by = PC7, bs = "cc") +
#                                  s(Doy, by = PC8, bs = "cc") +
#                                  s(Doy, by = PC9, bs = "cc") +
#                                  s(Doy, by = PC10, bs = "cc") +
#                                  s(Doy, by = PC11, bs = "cc") +
#                                  s(Doy, by = PC12, bs = "cc"), # +
#                                # s(Doy, by = PC13, bs = "cc") +
#                                # s(Doy, by = PC14, bs = "cc") +
#                                # s(Doy, by = PC15, bs = "cc") +
#                                # s(Doy, by = PC16, bs = "cc") +
#                                # s(Doy, by = PC17, bs = "cc") +
#                                # s(Doy, by = PC18, bs = "cc") +
#                                # s(Doy, by = PC19, bs = "cc") +
#                                # s(Doy, by = PC20, bs = "cc") +
#                                # s(Doy, by = PC21, bs = "cc") +
#                                # s(Doy, by = PC22, bs = "cc"),
#                                data = all_dat[train_id,])
# # plot(lm_ab_annualcycle2, scheme = 2, se = F)
# anova(lm_ab_annualcycle2)
# 
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# cor(preds,
#     all_dat[test_id, "King_Hakon"], use = "pairwise")
# rmse = function(x1, x2) {
#   sqrt(mean((x1 - x2)^2, na.rm = TRUE))
# }
# plot(all_dat$Date[test_id], 
#      all_dat[test_id, "King_Hakon"], type = "l")
# lines(all_dat$Date[test_id], preds, col = 2)
# preds2 = predict(lm_ab_annualcycle2, 
#                  all_dat[test_id,])
# cor(preds2,
#     all_dat[test_id, "King_Hakon"], use = "pairwise")
# lines(all_dat$Date[test_id], preds2, col = 3)
# abline(h = 0, lty = 2)
# rmse(0, all_dat[test_id, "King_Hakon"])
# rmse(preds, all_dat[test_id, "King_Hakon"])
# rmse(preds2, all_dat[test_id, "King_Hakon"])
# 
# # Ross
# lm_ab_annualcycle = mgcv::gam(Ross ~ te(Doy, PC1, bs = c("cc", "ds")) +
#                                 te(Doy, PC2, bs = c("cc", "ds")) +
#                                 te(Doy, PC3, bs = c("cc", "ds")) +
#                                 te(Doy, PC4, bs = c("cc", "ds")) +
#                                 te(Doy, PC5, bs = c("cc", "ds")), 
#                               data = all_dat[train_id,])
# # plot(lm_ab_annualcycle, scheme = 2, se = F)
# lm_ab_annualcycle2 = mgcv::gam(Ross ~ s(Doy, by = PC1, bs = "cc") +
#                                  s(Doy, by = PC2, bs = "cc") +
#                                  s(Doy, by = PC3, bs = "cc") +
#                                  s(Doy, by = PC4, bs = "cc") +
#                                  s(Doy, by = PC5, bs = "cc") +
#                                  s(Doy, by = PC6, bs = "cc") +
#                                  s(Doy, by = PC7, bs = "cc") +
#                                  s(Doy, by = PC8, bs = "cc") +
#                                  s(Doy, by = PC9, bs = "cc") +
#                                  s(Doy, by = PC10, bs = "cc") +
#                                  s(Doy, by = PC11, bs = "cc") +
#                                  s(Doy, by = PC12, bs = "cc"), #  +
#                                # s(Doy, by = PC13, bs = "cc") +
#                                # s(Doy, by = PC14, bs = "cc") +
#                                # s(Doy, by = PC15, bs = "cc") +
#                                # s(Doy, by = PC16, bs = "cc") +
#                                # s(Doy, by = PC17, bs = "cc") +
#                                # s(Doy, by = PC18, bs = "cc") +
#                                # s(Doy, by = PC19, bs = "cc") +
#                                # s(Doy, by = PC20, bs = "cc") +
#                                # s(Doy, by = PC21, bs = "cc") +
#                                # s(Doy, by = PC22, bs = "cc"), 
#                                data = all_dat[train_id,])
# # plot(lm_ab_annualcycle2, scheme = 2, se = F)
# anova(lm_ab_annualcycle2)
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# cor(preds,
#     all_dat[test_id, "Ross"], use = "pairwise")
# preds2 = predict(lm_ab_annualcycle2, 
#                  all_dat[test_id,])
# cor(preds2,
#     all_dat[test_id, "Ross"], use = "pairwise")
# plot(all_dat$Date[test_id], 
#      all_dat[test_id, "Ross"], type = "l")
# lines(all_dat$Date[test_id], preds, col = 2)
# lines(all_dat$Date[test_id], preds2, col = 3)
# abline(h = 0, lty = 2)
# rmse(0, all_dat[test_id, "Ross"])
# rmse(preds, all_dat[test_id, "Ross"])
# rmse(preds2, all_dat[test_id, "Ross"])
# 
# lm_ab_annualcycle = mgcv::gam(East_Antarctica ~ te(Doy, PC1, bs = c("cc", "ps")) +
#                                 te(Doy, PC2, bs = c("cc", "ds")) +
#                                 te(Doy, PC3, bs = c("cc", "ds")) +
#                                 te(Doy, PC4, bs = c("cc", "ds")) +
#                                 te(Doy, PC5, bs = c("cc", "ds")), 
#                               data = all_dat[train_id,])
# # plot(lm_ab_annualcycle, scheme = 2, se = F)
# lm_ab_annualcycle2 = mgcv::gam(East_Antarctica ~ s(Doy, by = PC1, bs = "cc") +
#                                  s(Doy, by = PC2, bs = "cc") +
#                                  s(Doy, by = PC3, bs = "cc") +
#                                  s(Doy, by = PC4, bs = "cc") +
#                                  s(Doy, by = PC5, bs = "cc") +
#                                  s(Doy, by = PC6, bs = "cc") +
#                                  s(Doy, by = PC7, bs = "cc") +
#                                  s(Doy, by = PC8, bs = "cc") +
#                                  s(Doy, by = PC9, bs = "cc") +
#                                  s(Doy, by = PC10, bs = "cc") +
#                                  s(Doy, by = PC11, bs = "cc") +
#                                  s(Doy, by = PC12, bs = "cc"), # +
#                                # s(Doy, by = PC13, bs = "cc") +
#                                # s(Doy, by = PC14, bs = "cc") +
#                                # s(Doy, by = PC15, bs = "cc") +
#                                # s(Doy, by = PC16, bs = "cc") +
#                                # s(Doy, by = PC17, bs = "cc") +
#                                # s(Doy, by = PC18, bs = "cc") +
#                                # s(Doy, by = PC19, bs = "cc") +
#                                # s(Doy, by = PC20, bs = "cc") +
#                                # s(Doy, by = PC21, bs = "cc") +
#                                # s(Doy, by = PC22, bs = "cc"),
#                                data = all_dat[train_id,])
# # plot(lm_ab_annualcycle2, scheme = 2, se = F)
# 
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# cor(preds,
#     all_dat[test_id, "East_Antarctica"], use = "pairwise")
# preds2 = predict(lm_ab_annualcycle2, 
#                  all_dat[test_id,])
# cor(preds2,
#     all_dat[test_id, "East_Antarctica"], use = "pairwise")
# 
# plot(all_dat$Date[test_id], 
#      all_dat[test_id, "East_Antarctica"], type = "l")
# lines(all_dat$Date[test_id], preds, col = 2)
# lines(all_dat$Date[test_id], preds2, col = 3)
# abline(h = 0, lty = 2)
# rmse(0, all_dat[test_id, "East_Antarctica"])
# rmse(preds, all_dat[test_id, "East_Antarctica"])
# rmse(preds2, all_dat[test_id, "East_Antarctica"])
# 
# lm_ab_annualcycle = mgcv::gam(Weddell ~ te(Doy, PC1, bs = c("cc", "ps")) +
#                                 te(Doy, PC2, bs = c("cc", "ds")) +
#                                 te(Doy, PC3, bs = c("cc", "ds")) +
#                                 te(Doy, PC4, bs = c("cc", "ds")) +
#                                 te(Doy, PC5, bs = c("cc", "ds")),
#                               data = all_dat[train_id,])
# # plot(lm_ab_annualcycle, scheme = 2, se = F)
# lm_ab_annualcycle2 = mgcv::gam(Weddell ~  s(Doy, by = PC1, bs = "cc") +
#                                  s(Doy, by = PC2, bs = "cc") +
#                                  s(Doy, by = PC3, bs = "cc") +
#                                  s(Doy, by = PC4, bs = "cc") +
#                                  s(Doy, by = PC5, bs = "cc") +
#                                  s(Doy, by = PC6, bs = "cc") +
#                                  s(Doy, by = PC7, bs = "cc") +
#                                  s(Doy, by = PC8, bs = "cc") +
#                                  s(Doy, by = PC9, bs = "cc") +
#                                  s(Doy, by = PC10, bs = "cc") +
#                                  s(Doy, by = PC11, bs = "cc") +
#                                  s(Doy, by = PC12, bs = "cc"), # +
#                                # s(Doy, by = PC13, bs = "cc") +
#                                # s(Doy, by = PC14, bs = "cc") +
#                                # s(Doy, by = PC15, bs = "cc") +
#                                # s(Doy, by = PC16, bs = "cc") +
#                                # s(Doy, by = PC17, bs = "cc") +
#                                # s(Doy, by = PC18, bs = "cc") +
#                                # s(Doy, by = PC19, bs = "cc") +
#                                # s(Doy, by = PC20, bs = "cc") +
#                                # s(Doy, by = PC21, bs = "cc") +
#                                # s(Doy, by = PC22, bs = "cc"), 
#                                data = all_dat[train_id,])
# # plot(lm_ab_annualcycle2, scheme = 2, se = F)
# 
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# cor(preds,
#     all_dat[test_id, "Weddell"], use = "pairwise")
# preds2 = predict(lm_ab_annualcycle2, 
#                  all_dat[test_id,])
# cor(preds2,
#     all_dat[test_id, "Weddell"], use = "pairwise")
# plot(all_dat$Date[test_id], 
#      all_dat[test_id, "Weddell"], type = "l")
# lines(all_dat$Date[test_id], preds, col = 2)
# lines(all_dat$Date[test_id], preds2, col = 3)
# abline(h = 0, lty = 2)
# rmse(0, all_dat[test_id, "Weddell"])
# rmse(preds, all_dat[test_id, "Weddell"])
# rmse(preds2, all_dat[test_id, "Weddell"])
# 
# lm_ab_annualcycle = mgcv::gam(Bellingshausen_Amundsen_Sea ~ te(Doy, PC1, bs = c("cc", "ps")) +
#                                 te(Doy, PC2, bs = c("cc", "ds")) +
#                                 te(Doy, PC3, bs = c("cc", "ds")) +
#                                 te(Doy, PC4, bs = c("cc", "ds")) +
#                                 te(Doy, PC5, bs = c("cc", "ds")), 
#                               data = all_dat[train_id,])
# # plot(lm_ab_annualcycle, scheme = 2, se = F)
# lm_ab_annualcycle2 = mgcv::gam(Bellingshausen_Amundsen_Sea ~ s(Doy, by = PC1, bs = "cc") +
#                                  s(Doy, by = PC2, bs = "cc") +
#                                  s(Doy, by = PC3, bs = "cc") +
#                                  s(Doy, by = PC4, bs = "cc") +
#                                  s(Doy, by = PC5, bs = "cc") +
#                                  s(Doy, by = PC6, bs = "cc") +
#                                  s(Doy, by = PC7, bs = "cc") +
#                                  s(Doy, by = PC8, bs = "cc") +
#                                  s(Doy, by = PC9, bs = "cc") +
#                                  s(Doy, by = PC10, bs = "cc") +
#                                  s(Doy, by = PC11, bs = "cc") +
#                                  s(Doy, by = PC12, bs = "cc") +
#                                s(Doy, by = PC13, bs = "cc") +
#                                s(Doy, by = PC14, bs = "cc") +
#                                s(Doy, by = PC15, bs = "cc") +
#                                s(Doy, by = PC16, bs = "cc") +
#                                s(Doy, by = PC17, bs = "cc") +
#                                s(Doy, by = PC18, bs = "cc") +
#                                s(Doy, by = PC19, bs = "cc") +
#                                s(Doy, by = PC20, bs = "cc"),
#                                # s(Doy, by = PC21, bs = "cc") +
#                                # s(Doy, by = PC22, bs = "cc"), 
#                                data = all_dat[train_id,])
# # plot(lm_ab_annualcycle2, scheme = 2, se = F)
# 
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# cor(preds,
#     all_dat[test_id, "Bellingshausen_Amundsen_Sea"], use = "pairwise")
# preds2 = predict(lm_ab_annualcycle2, 
#                  all_dat[test_id,])
# cor(preds2,
#     all_dat[test_id, "Bellingshausen_Amundsen_Sea"], use = "pairwise")
# plot(all_dat$Date[test_id], 
#      all_dat[test_id, "Bellingshausen_Amundsen_Sea"], type = "l")
# lines(all_dat$Date[test_id], preds, col = 2)
# lines(all_dat$Date[test_id], preds2, col = 3)
# abline(h = 0, lty = 2)
# rmse(0, all_dat[test_id, "Bellingshausen_Amundsen_Sea"])
# rmse(preds, all_dat[test_id, "Bellingshausen_Amundsen_Sea"])
# rmse(preds2, all_dat[test_id, "Bellingshausen_Amundsen_Sea"])


# # ########################################
# # combine lags and nonlinear relationships
# lags = c(2, 7, 29, 61)
# i = 1
# for (i in seq_along(lags)) {
#   all_dat[,paste0("PC1_l", lags[i])] = lag(all_dat$PC1, n = lags[i])
#   all_dat[,paste0("PC2_l", lags[i])] = lag(all_dat$PC1, n = lags[i])
#   all_dat[,paste0("PC3_l", lags[i])] = lag(all_dat$PC1, n = lags[i])
#   all_dat[,paste0("PC4_l", lags[i])] = lag(all_dat$PC1, n = lags[i])
#   all_dat[,paste0("PC5_l", lags[i])] = lag(all_dat$PC1, n = lags[i])
# }
# 
# res = data.frame(lag = lags)
# # for (i in seq_along(lags)) {
# # TODO there are NA is long lags due to data 
# lm_ab_annualcycle = mgcv::gam(King_Hakon ~ te(Doy, PC1, bs = c("cc", "ds")) +
#                                 te(Doy, PC1_l2, bs = c("cc", "ds")) +
#                                 te(Doy, PC1_l7, bs = c("cc", "ds")) +
#                                 # te(Doy, PC1_l17, bs = c("cc", "ds")) +
#                                 te(Doy, PC1_l29, bs = c("cc", "ds")) +
#                                 te(Doy, PC1_l61, bs = c("cc", "ds")) +
#                                 te(Doy, PC2, bs = c("cc", "ds")) +
#                                 # te(Doy, PC2_l7, bs = c("cc", "ds")) +
#                                 # te(Doy, PC2_l30, bs = c("cc", "ds")) +
#                                 te(Doy, PC3, bs = c("cc", "ds")) +
#                                 # te(Doy, PC3_l7, bs = c("cc", "ds")) +
#                                 te(Doy, PC4, bs = c("cc", "ds")) +
#                                 # te(Doy, PC4_l7, bs = c("cc", "ds")) +
#                                 te(Doy, PC5, bs = c("cc", "ds")), # +
#                               # te(Doy, PC5_l7, bs = c("cc", "ds")), 
#                               data = all_dat)
# anova(lm_ab_annualcycle)
# plot(lm_ab_annualcycle, scheme = 2, se = F)
# # plot(test_id, all_dat$PC1_lagged[test_id], type = "l")
# preds = predict(lm_ab_annualcycle,
#                 all_dat[test_id,])
# res$King_Hakon_oob[i] = mean((preds - all_dat[test_id, "King_Hakon"])^2, na.rm = TRUE)
# res$King_Hakon_sse[i] = mean(lm_ab_annualcycle$residuals^2)
# cor(preds,
#     all_dat[test_id, "King_Hakon"], use = "pairwise")
# plot(test_id,
#      all_dat[test_id, "King_Hakon"], type = "l")
# lines(preds, col = 2)
# 
# # Ross
# lm_ab_annualcycle = mgcv::gam(Ross ~ te(Doy, PC1_lagged, bs = c("cc", "ds")), data = all_dat)
# plot(lm_ab_annualcycle, scheme = 2, se = F)
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# res$Ross_oob[i] = mean((preds - all_dat[test_id, "Ross"])^2, na.rm = TRUE)
# res$Ross_sse[i] = mean(lm_ab_annualcycle$residuals^2)
# # # why the negative correlation?
# # cor(preds,
# #     all_dat[test_id, "Ross"], use = "pairwise")
# # plot(test_id,
# #      all_dat[test_id, "Ross"], type = "l")
# # lines(preds, col = 2)
# # abline(h = 0, lty = 2)
# 
# lm_ab_annualcycle = mgcv::gam(East_Antarctica ~ te(Doy, PC1_lagged, bs = c("cc", "ps")), data = all_dat)
# plot(lm_ab_annualcycle, scheme = 2, se = F)
# # lm_ab_annualcycle
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# res$East_Antarctica_oob[i] = mean((preds - all_dat[test_id, "Ross"])^2, na.rm = TRUE)
# res$East_Antarctica_sse[i] = mean(lm_ab_annualcycle$residuals^2)
# # cor(preds,
# #     all_dat[test_id, "East_Antarctica"], use = "pairwise")
# # plot(test_id,
# #      all_dat[test_id, "East_Antarctica"], type = "l")
# # lines(preds, col = 2)
# # abline(h = 0, lty = 2)
# 
# lm_ab_annualcycle = mgcv::gam(Weddell ~ te(Doy, PC1_lagged, bs = c("cc", "ps")), data = all_dat)
# plot(lm_ab_annualcycle, scheme = 2, se = F)
# # lm_ab_annualcycle
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# res$Weddell_oob[i] = mean((preds - all_dat[test_id, "Weddell"])^2, na.rm = TRUE)
# res$Weddell_sse[i] = mean(lm_ab_annualcycle$residuals^2)
# # cor(preds,
# #     all_dat[test_id, "Weddell"], use = "pairwise")
# # plot(test_id,
# #      all_dat[test_id, "Weddell"], type = "l")
# # lines(preds, col = 2)
# # abline(h = 0, lty = 2)
# 
# lm_ab_annualcycle = mgcv::gam(Bellingshausen_Amundsen_Sea ~ te(Doy, PC1_lagged, bs = c("cc", "ps")), data = all_dat)
# plot(lm_ab_annualcycle, scheme = 2, se = F)
# # lm_ab_annualcycle
# preds = predict(lm_ab_annualcycle, 
#                 all_dat[test_id,])
# res$Bellingshausen_Amundsen_Sea_oob[i] = mean((preds - all_dat[test_id, "Weddell"])^2, na.rm = TRUE)
# res$Bellingshausen_Amundsen_Sea_sse[i] = mean(lm_ab_annualcycle$residuals^2)
# # cor(preds,
# #     all_dat[test_id, "Bellingshausen_Amundsen_Sea"], use = "pairwise")
# # plot(test_id,
# #      all_dat[test_id, "Bellingshausen_Amundsen_Sea"], type = "l")
# # lines(preds, col = 2)
# # abline(h = 0, lty = 2)
# # }
# # in bag predictions
# plot(lags, res$King_Hakon_sse, type = "l", ylim = c(0.1, 0.14))
# lines(lags, res$Ross_sse)
# lines(lags, res$Weddell_sse)
# plot(lags, res$East_Antarctica_sse, type = "l", ylim = c(0.03, 0.036))
# lines(lags, res$Bellingshausen_Amundsen_Sea_sse)
# 
# # out of bag predictions
# plot(lags, res$King_Hakon_oob, type = "l", ylim = c(0.0, 0.2))
# lines(lags, res$Ross_oob)
# lines(lags, res$Weddell_oob)
# 
# plot(lags, res$East_Antarctica_oob, type = "l", ylim = c(0.0, 0.2))
# lines(lags, res$Bellingshausen_Amundsen_Sea_oob)



# # create linear model on PC level
# target_amp = t(t(eigen_Y$rotation) %*% Y)
# cov_amp = t(t(V) %*% X)
# 
# lags = seq(0, 360, by = 5)
# r2_pc1 = rep(NA, length(lags))
# r2_pc2 = rep(NA, length(lags))
# r2_pc3 = rep(NA, length(lags))
# r2_pc4 = rep(NA, length(lags))
# r2_pc5 = rep(NA, length(lags))
# 
# for (i in seq_along(lags)) {
#   cov_amp_lag = apply(cov_amp, 2, lag, n = lags[i])
#   lm_ab = lm(target_amp ~ cov_amp_lag)
#   lm_ab_sum = summary(lm_ab)
#   r2_pc1[i] = lm_ab_sum$`Response PC1`$r.squared
#   r2_pc2[i] = lm_ab_sum$`Response PC2`$r.squared
#   r2_pc3[i] = lm_ab_sum$`Response PC3`$r.squared
#   r2_pc4[i] = lm_ab_sum$`Response PC4`$r.squared
#   r2_pc5[i] = lm_ab_sum$`Response PC5`$r.squared
# }
# 
# plot(lags, r2_pc1, ylim = c(0, 0.05), type = "l", ylab = "r squared")
# lines(lags, r2_pc2, col = 2)
# lines(lags, r2_pc3, col = 3)
# lines(lags, r2_pc4, col = 3)
# lines(lags, r2_pc5, col = 3)
# 
# 
# # create linear model on actual level
# target_amp = t(Y)
# cov_amp = t(t(V) %*% X)
# 
# lags = seq(0, 360, by = 5)
# r2_1 = rep(NA, length(lags))
# r2_2 = rep(NA, length(lags))
# r2_3 = rep(NA, length(lags))
# r2_4 = rep(NA, length(lags))
# r2_5 = rep(NA, length(lags))

# # TODO run separately for four seasons or months
# # or plug in mgcv::gam from below
# # or try sparse regression (LASSO) on all lags at the same time
# for (i in seq_along(lags)) {
#   cov_amp_lag = apply(cov_amp, 2, lag, n = lags[i])
#   lm_ab = lm(target_amp[c(rep(F, 120), rep(TRUE, 30), rep(F, 215)),] ~ 
#                cov_amp_lag[c(rep(F, 120), rep(TRUE, 30), rep(F, 215)),])
#   lm_ab_sum = summary(lm_ab)
#   r2_1[i] = lm_ab_sum$`Response King_Hakon`$r.squared
#   r2_2[i] = lm_ab_sum$`Response Ross`$r.squared
#   r2_3[i] = lm_ab_sum$`Response East_Antarctica`$r.squared
#   r2_4[i] = lm_ab_sum$`Response Weddell`$r.squared
#   r2_5[i] = lm_ab_sum$`Response Bellingshausen_Amundsen_Sea`$r.squared
# }
# 
# plot(lags, r2_1, ylim = c(0, 0.3), type = "l", ylab = "r squared")
# lines(lags, r2_2, col = 2)
# lines(lags, r2_3, col = 3)
# lines(lags, r2_4, col = 4)
# lines(lags, r2_5, col = 5)


