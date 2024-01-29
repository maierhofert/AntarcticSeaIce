# manually do something like a PCA
# just aggregate the variables geospatially
library("dplyr")

# analyze the combined covariates available
source("cov_daily.R")
source("cov_downsample_monthly_to_daily.R")
# TODO there is a SST reconstruction that I am not using currently
# source("cov_climate_indices.R")

# combine all covariates
cov = full_join(daily_cov_daily, mon_cov_daily, by = "Date")
# cov = full_join(cov, ind_cov_daily, by = "Date")
cov = cov[order(cov$Date),]
summary(is.na(cov))

# 2019 is not available for monthly data, but for daily and for indices
# we also do not have the sea ice in 2019, so that should be alright for now
cov$Date[which(is.na(cov$'619010_SLP'))]
cov = cov[cov$Date < as.Date("2019-01-01"),]

names(cov)

# check out sd
sds = apply(cov[,-1], 2, sd)
hist(sds)
# check out means
means = colMeans(cov[,-1])
hist(means)

# standardize
cov2 = cov
for (i in 2:ncol(cov)) {
  cov2[,i] = (cov[,i] - means[i-1]) / sds[i-1]
}
cov = cov2


# get groups for aggrgeation
names(cov)

groupings = list(TMINMAX_0 = c("TMAX_ASN00002011", "TMIN_ASN00002011", "TMAX_ASN00003002", "TMIN_ASN00003002",
                               "TMAX_ASN00004020", "TMIN_ASN00004020",
                               "TMAX_ASN00009034", "TMIN_ASN00009034", "TMAX_ASN00010648", "TMIN_ASN00010648"),
                 TMINMAX_1_2 = c("TMAX_ASN00015540", "TMIN_ASN00015540",
                                 "TMAX_ASN00018070", "TMIN_ASN00018070",
                                 "TMAX_ASN00023000", "TMIN_ASN00023000", "TMAX_ASN00026026", "TMIN_ASN00026026",
                                 "TMAX_ASN00028004", "TMIN_ASN00028004"), 
                 TMINMAX_3 = c("TMAX_ASN00030018", "TMIN_ASN00030018", "TAVG_ASN00030045", "TMAX_ASN00030045", "TMIN_ASN00030045",
                               "TAVG_ASN00031010", "TMAX_ASN00031010", "TMIN_ASN00031010", "TMAX_ASN00033001",
                               "TMIN_ASN00033001", "TMAX_ASN00033046", "TMIN_ASN00033046", "TMAX_ASN00034002", "TMIN_ASN00034002",
                               "TMAX_ASN00035027", "TMIN_ASN00035027", "TMAX_ASN00036030", "TMIN_ASN00036030", "TMAX_ASN00038003",
                               "TMIN_ASN00038003", "TMAX_ASN00039015", "TMIN_ASN00039015", "TMAX_ASN00039039", "TMIN_ASN00039039"),
                 TMINMAX_4_5_6 = c("TMAX_ASN00040214", "TMIN_ASN00040214", "TMAX_ASN00040264", "TMIN_ASN00040264", "TMAX_ASN00041023",
                                   "TMIN_ASN00041023", "TMAX_ASN00041038", "TMIN_ASN00041038", "TMAX_ASN00044022", "TMIN_ASN00044022",
                                   "TMAX_ASN00048013", "TMIN_ASN00048013", "TMAX_ASN00048030", "TMIN_ASN00048030", 
                                   "TMAX_ASN00052026", "TMIN_ASN00052026", "TMAX_ASN00055023", "TMIN_ASN00055023", "TMAX_ASN00056017", "TMIN_ASN00056017",    
                                   "TMAX_ASN00063004", "TMIN_ASN00063004", "TMAX_ASN00065016", "TMIN_ASN00065016", "TMAX_ASN00066062",      
                                   "TMIN_ASN00066062"), 
                 TMINMAX_7_8_9 = c("TMAX_ASN00072151", "TMIN_ASN00072151", "TMAX_ASN00074128", "TMIN_ASN00074128",
                                   "TMAX_ASN00076077", "TMIN_ASN00076077", "TMAX_ASN00078031", "TMIN_ASN00078031",
                                   "TMAX_ASN00086071", "TMIN_ASN00086071", "TAVG_ASN00090015", "TMAX_ASN00090015", "TMIN_ASN00090015", "TMAX_ASN00091049",      
                                   "TMIN_ASN00091049", "TMAX_ASN00091057", "TMIN_ASN00091057", "TMAX_ASN00094029", "TMIN_ASN00094029"),
                 #
                 TMP_6 = c("619010_tmp", "619960_tmp", "619980_tmp", "685880_tmp", "688160_tmp", "688420_tmp", "688580_tmp","689060_tmp",  "689940_tmp"), 
                 SLP_6 = c("619010_SLP", "619960_SLP", "619980_SLP", "685880_SLP", "688160_SLP", "688420_SLP", "688580_SLP", "689060_SLP","689940_SLP" ),
                 SLP_83_87 = c("837430_SLP", "854420_SLP", "854880_SLP", "855740_SLP", "855850_SLP", "857660_SLP",
                               "859340_SLP", "870470_SLP", "871550_SLP", "872220_SLP", 
                               "873440_SLP", "874800_SLP", "875440_SLP", "875480_SLP", "875760_SLP", "875850_SLP",
                               "876230_SLP", "876920_SLP", "877150_SLP", "877500_SLP", "878280_SLP", "878490_SLP",
                               "878600_SLP", "879830_SLP"),
                 TMP_83_87 = c("837430_tmp", "854420_tmp", "854880_tmp", "855740_tmp", "855850_tmp",  
                               "859340_tmp", "870470_tmp",  "871550_tmp", "872220_tmp", "872720_tmp", 
                               "873440_tmp", "874800_tmp", "875440_tmp",  "875480_tmp", "875760_tmp",  "875850_tmp",
                               "876230_tmp", "876920_tmp", "877150_tmp", "877500_tmp", "878280_tmp",  "878490_tmp", 
                               "878600_tmp", "879830_tmp"), 
                 SLP_88_89 = c("889030_SLP", "889630_SLP", "889680_SLP", "890090_SLP", "890220_SLP",
                               "890500_SLP","890550_SLP", "890620_SLP", "890630_SLP", "891250_SLP",
                               "895120_SLP", "895320_SLP", "895640_SLP", "895710_SLP",
                               "895920_SLP", "896060_SLP", "896110_SLP", "896420_SLP"),
                 TMP_88_89 = c("889030_tmp", "889630_tmp", "889680_tmp", "890090_tmp", "890220_tmp", 
                               "890500_tmp", "890550_tmp", "890620_tmp", "890630_tmp", 
                               "895120_tmp", "895320_tmp", "895640_tmp", "895710_tmp",  
                               "895920_tmp", "896060_tmp", "896110_tmp", "896420_tmp"),
                 "Marsh_O_Higgins_SLP", "Marsh_O_Higgins_tmp", "Scott_Base_mcMurdo_SLP", "Scott_Base_McMurdo_tmp", 
                 SLP_91_93 = c("915770_SLP", "915920_SLP", "916800_SLP",  "918430_SLP",
                               "919380_SLP", "931190_SLP", "934340_SLP", "936150_SLP",
                               "937800_SLP", "938940_SLP", "939470_SLP", "939870_SLP"),
                 TMP_91_93 = c("915770_tmp",  "915920_tmp", "916430_tmp", "916800_tmp", "918430_tmp",
                               "919380_tmp", "931190_tmp", "934340_tmp", "936150_tmp", 
                               "937800_tmp", "938940_tmp", "939470_tmp", "939870_tmp"), 
                 SLP_94 = c("943000_SLP", "943120_SLP", "943260_SLP", "943460_SLP",
                            "943670_SLP", "944300_SLP", "945100_SLP", "945780_SLP", "946100_SLP",
                            "946370_SLP", "946530_SLP", "946720_SLP", "947680_SLP", "948680_SLP",
                            "949260_SLP", "949750_SLP", "949980_SLP"),
                 TMP_94 = c("943000_tmp", "943120_tmp", "943260_tmp", "943460_tmp",
                            "943670_tmp", "944300_tmp", "945100_tmp", "945780_tmp", "946100_tmp",    
                            "946370_tmp", "946530_tmp", "946720_tmp", "947680_tmp", "948680_tmp", 
                            "949260_tmp", "949750_tmp", "949960_tmp", "949980_tmp"))
# check correct number of variables
length(unlist(groupings))
length(unique(unlist(groupings)))
dim(cov)

cov_red = data.frame(Date = cov$Date)
for (i in 1:length(groupings)) {
  if (length(groupings[[i]]) == 1) {
    cov_red[,groupings[[i]]] = cov[,groupings[[i]]]
  } else {
    cov_red[,names(groupings[i])] = apply(cov[,groupings[[i]]], 1, mean) 
  }
}
head(cov_red)
dim(cov_red)

# compute correlation matrix
cor_cov = cor(cov_red[,-1])
# plot correlation matrix
library(ggplot2)
melted_cor_cov <- reshape2::melt(cor_cov)
ggplot(data = melted_cor_cov, aes(x=Var1, y=Var2, fill=value)) + 
  xlab("") + ylab("") +
  scale_fill_gradient2("Correlation", limits = c(-1, 1), high = "firebrick", low = "steelblue") +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../plots/covariates/all_covariates/corr_manual_cluster.png", device = "png", width = 25, height = 25)


for (i in 2:ncol(cov_red)) {
  plot(cov_red$Date[cov_red$Date > as.Date("2000-01-01")],
       cov_red[,i][cov_red$Date > as.Date("2000-01-01")],
       main = names(cov_red)[i], type = "l")
}
cov_red$Doy = lubridate::yday(cov_red$Date)

# test = mgcv::gam(TMINMAX_0 ~ s(Doy, bs = 'cc'), 
#           data = cov_red[cov_red$Date > as.Date("1900-01-01"),])
# plot(test)
# test = mgcv::gam(TMP_6 ~ s(Doy, bs = 'cc'), 
#                  data = cov_red)
# plot(test)
# test = mgcv::gam(SLP_6 ~ s(Doy, bs = 'cc'), 
#                  data = cov_red)
# plot(test)
# test = mgcv::gam(Marsh_O_Higgins_SLP ~ s(Doy, bs = 'cc'), 
#                  data = cov_red[cov_red$Date > as.Date("1950-01-01"),])
# plot(test)
# test = mgcv::gam(Marsh_O_Higgins_tmp ~ s(Doy, bs = 'cc'), 
#                  data = cov_red)
# plot(test)

# standardize the variables
sds = apply(cov_red, 2, sd)
sds
means = colMeans(cov_red[,-1])
means %>% round(2)

for (i in 2:(ncol(cov_red) - 1)) {
  cov_red[,i] = cov_red[,i] / sds[i]
}
