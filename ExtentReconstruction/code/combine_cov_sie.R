# combine all data
# the following files create the data sets
# read_Fogt_predictor_monthly.R impute_monthly_amelia.R cov_monthly.R target_sea_ice.R

# read covariates
cov = readRDS("../data/fogt_predictor_data_imputed/all_monthly_cov.RDS")
covLag = readRDS("../data/fogt_predictor_data_imputed/all_monthly_cov_lagged.RDS")
covPC = readRDS("../data/fogt_predictor_data_imputed/all_monthly_PC_lagged.RDS")

# read target variable
sie = readRDS("../data/HandockRegionalSeaIce/sie_detrended_mon.RDS")
sie$Date = as.Date(paste0(sie$Year, "-", sie$Month, "-15"))

# combine the covariates and target sea ice
all_data = full_join(sie, cov, by = c("Year", "Month", "Date"))
all_data = all_data[order(all_data$Date),]
saveRDS(all_data, "../data/complete_dataset/all_data.RDS")

all_dataLag = full_join(sie, covLag, by = c("Year", "Month", "Date"))
all_dataLag = all_dataLag[order(all_dataLag$Date),]
saveRDS(all_dataLag, "../data/complete_dataset/all_dataLag.RDS")

all_dataPC = full_join(sie, covPC, by = c("Year", "Month", "Date"))
all_dataPC = all_dataPC[order(all_dataPC$Date),]
saveRDS(all_dataPC, "../data/complete_dataset/all_dataPC.RDS")

