library("dplyr")
# load sea ice data
load("../data/HandockRegionalSeaIce/nsidcV4_Regions.RData")

# load("../data/HandockRegionalSeaIce/nsidcG02202_V4_Regions.RData")

month.m = round((tdate.m %% 1) * 12 + 0.5)
colnames(sie.m)[5] = "Bellingshausen_Amundsen_Sea"
all_data = data.frame(tdate.m, Month = month.m, Year = year.m, sie.m)

# looking good
plot(tdate.m, sie.m[,1], type = "l")
plot(tdate.m, sie.m[,2], type = "l")
plot(tdate.m, sie.m[,3], type = "l")
plot(tdate.m, sie.m[,4], type = "l")
plot(tdate.m, sie.m[,5], type = "l")


# total sea ice
load("../data/HandockRegionalSeaIce/nsidcV4.RData")

load("../data/HandockRegionalSeaIce/nsidcV42023.RData")
my_dat = data.frame(sie.m = sie.m, tdate.m = tdate.m, Month = tdate.m %% 1)
res = my_dat %>% group_by(Month) %>%
  summarize(Extent = mean(sie.m, na.rm = TRUE))
res
tail(my_dat, 12)

# daily data
my_dat.daily = data.frame(sie = sie, tdate = tdate, day = tdate %% 1)
tail(my_dat.daily)
Feb2023 = my_dat.daily[my_dat.daily$tdate > 2023.086 & my_dat.daily$tdate < 2023.161,]
dim(Feb2023)
mean(Feb2023$sie)

# load("../data/HandockRegionalSeaIce/nsidcG02202_V4.RData")
# sie.m = sie.m[-(529:540)]
all_data = cbind(all_data, total = sie.m)
all_data %>% group_by(Month) %>%
  summarize(Extent = mean(total, na.rm = TRUE))
# plot(tdate.m, sie.m)

# check sea ice
plot(all_data$tdate.m, all_data$King_Hakon)
plot(all_data$tdate.m, all_data$Ross)
plot(all_data$tdate.m, all_data$East_Antarctica)
plot(all_data$tdate.m, all_data$Weddell)
plot(all_data$tdate.m, all_data$Bellingshausen_Amundsen_Sea)
plot(all_data$tdate.m, all_data$total)

# Decmber 1987 was wrong, it's fixed now
# pdf("../plots/SeaIce1987.pdf", height = 4, width = 7)
plot(all_data$tdate.m, all_data$King_Hakon, xlim = c(1987,1989))
plot(all_data$tdate.m, all_data$Ross, xlim = c(1987,1989))
plot(all_data$tdate.m, all_data$East_Antarctica, xlim = c(1987,1989))
plot(all_data$tdate.m, all_data$Weddell, xlim = c(1987,1989))
plot(all_data$tdate.m, all_data$Bellingshausen_Amundsen_Sea, xlim = c(1987,1989))
plot(all_data$tdate.m, all_data$total, xlim = c(1987,1989))
dev.off()


# the numbers match except for a few months with missing values, but that's all good
plot(rowSums(all_data[,4:8]), all_data$total)
plot(all_data$tdate, rowSums(all_data[,4:8]) -  all_data$total)

sie = all_data[,4:8]
eigen_sie = prcomp(sie[complete.cases(sie),], scale. = TRUE, center = TRUE)
plot(eigen_sie$sdev / sum(eigen_sie$sdev))

mod_data = all_data[,c(2, 4:9)]
sie_detrended = mod_data

# detrend sea ice
for (i in 2:7) {
  this.var = colnames(mod_data)[i]
  this.gam = mgcv::gam(formula(paste0(this.var, " ~ s(Month, bs = 'cc')")),
                       knots = list(Month = c(0.5, 12.5)),
                       data = mod_data)
  sie_detrended[,this.var] = mod_data[,this.var] - predict(this.gam, newdata = mod_data)
  # TODO should scale be TRUE or FALSE
  # Do we care about the different scale of the different sectors
  # sie_detrended[,this.var] = scale(sie_detrended[,this.var], center = TRUE, scale = FALSE)
}
sie_detrended$Year = year.m
sie_detrended = sie_detrended[,c("Year", "Month", "total", "King_Hakon", "Ross",
                                 "East_Antarctica", "Weddell", "Bellingshausen_Amundsen_Sea")]
saveRDS(sie_detrended, "../data/HandockRegionalSeaIce/sie_detrended_mon.RDS")

sie_detrended = readRDS("../data/HandockRegionalSeaIce/sie_detrended_mon.RDS")
# February 2022 is 2.224553 with an anomaly of -0.9086785543
# February avg 1979-2022 is 3.133231

# February 2023 is 2.044898 with an anomaly of 1.088333
# February avg 1979-2022 is 3.133231
plot(sie_detrended$Year + (sie_detrended$Month - 0.5) / 12,
     sie_detrended$total)
plot(sie_detrended$Year + (sie_detrended$Month - 0.5) / 12,
     sie_detrended$King_Hakon)

all_data = all_data[,c("Year", "Month", "total", "King_Hakon", "Ross",
                       "East_Antarctica", "Weddell", "Bellingshausen_Amundsen_Sea")]
saveRDS(all_data, "../data/HandockRegionalSeaIce/sie_mon.RDS")


# MH use an arcsin trafo to standardize?
eigen_sie_detrended = prcomp(sie_detrended[complete.cases(sie_detrended),4:8], scale. = FALSE, center = TRUE)
plot(eigen_sie_detrended$sdev^2 / sum(eigen_sie_detrended$sdev^2))
eigen_sie_detrended$rotation %>% round(2)

# compare original variables var to PC
(eigen_sie_detrended$sdev^2 / sum(eigen_sie_detrended$sdev^2)) %>% round(2)
(apply(sie_detrended[,4:8], 2, sd, na.rm = TRUE) / sum(apply(sie_detrended[,4:8], 2, sd, na.rm = TRUE))) %>% round(2)

plot(sie_detrended[,3], type = "l")
plot(sie_detrended[,4], type = "l")
plot(sie_detrended[,5], type = "l")
plot(sie_detrended[,6], type = "l")
plot(sie_detrended[,7], type = "l")
plot(sie_detrended[,8], type = "l")

sie_date = as.Date(paste0(year, "-", month, "-", day))


# ######################################
# do an eigenvector decomposition on target

# compute correlation matrix
sie_detrended2 = sie_detrended
colnames(sie_detrended2) = c("Year", "Month", "Total", "King Hakon VII", "Ross", "East Antarctica",
                             "Weddell", "Bellingshausen\nAmundsen Sea")
cor_tar = cor(sie_detrended2, use = "pairwise")
# cor_tar = cov(sie_detrended2, use = "pairwise")
# plot correlation matrix
library(ggplot2)
melted_cor_tar <- reshape2::melt(cor_tar)
ggplot(data = melted_cor_tar, aes(x=Var1, y=Var2, fill=value)) +
  xlab("") + ylab("") +
  scale_fill_gradient2("Correlation", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  geom_text(aes(x = Var1, y = Var2, label = round(value, 2)))
ggsave("../plots/sea_ice/corr.png", device = "png", width = 7, height = 4)
# ggplot(data = melted_cor_tar, aes(x=Var1, y=Var2, fill=abs(value))) +
#   xlab("") + ylab("") +
#   scale_fill_gradient("abs(Correlation)", limits = c(0, 1),
#                       high = "firebrick", low = "white") +
#   geom_tile()
# ggsave("../plots/covariates/sea_ice/abs_corr.png", device = "png", width = 25, height = 25)

# plot the eigenvectors
eigen_sie_detrended = prcomp(sie_detrended2[complete.cases(sie_detrended),4:8], scale. = FALSE, center = TRUE)
melted_eigenvec <- reshape2::melt(eigen_sie_detrended$rotation)
ggplot(data = melted_eigenvec, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile() +
  xlab("") + ylab("") +
  scale_fill_gradient2("Loading", high = "firebrick", low = "steelblue", limits = c(-0.9, 0.9), oob = scales::squish) +
  geom_text(aes(x = Var2, y = Var1, label = round(value, 2)))
ggsave("../plots/sea_ice/eigenvec_col.png", device = "png", width = 7, height = 4)

ggplot(data = melted_eigenvec, aes(x=Var2, y=Var1, fill=abs(value))) +
  geom_tile() +
  xlab("") + ylab("") +
  scale_fill_gradient("abs(Loading)", high = "black", low = "white", limits = c(0, 0.9), oob = scales::squish) +
  geom_text(aes(x = Var2, y = Var1, label = round(value, 2)), color = "steelblue3")
ggsave("../plots/sea_ice/eigenvec_bw.png", device = "png", width = 7, height = 4)


