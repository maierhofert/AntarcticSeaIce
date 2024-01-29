my_DE_kriged_SLP = readRDS("../output/my_DE_kriged_SLP.RDS")
my_DE_kriged = readRDS("../output/my_DE_kriged.RDS")

# my_DE_kriged_SLP = readRDS("../output/my_DE_kriged_SLP_fineGrid.RDS")
# my_DE_kriged = readRDS("../output/my_DE_kriged_fineGrid.RDS")

test1 = data.frame(val = my_DE_kriged@data[, 1], 
                   variable = "TMP",
                   lon = my_DE_kriged@sp@coords[,1],
                   lat = my_DE_kriged@sp@coords[,2],
                   date = rep(spacetime:::trimDates(my_DE_kriged), 
                              each = nrow(my_DE_kriged@sp@coords)))

test2 = data.frame(val = my_DE_kriged_SLP@data[, 1], 
                   variable = "SLP",
                   lon = my_DE_kriged_SLP@sp@coords[,1],
                   lat = my_DE_kriged_SLP@sp@coords[,2],
                   date = rep(spacetime:::trimDates(my_DE_kriged_SLP), 
                              each = nrow(my_DE_kriged_SLP@sp@coords)))

kriged_data_long = rbind(test1, test2)
saveRDS(kriged_data_long, "../data/complete_dataset/gridded_long.RDS")

# ##############################################################################
# create wide data set
kriged_data_wide = data.frame(date = spacetime:::trimDates(my_DE_kriged))
# wide variables temperature at lon and lat per day
TMP_var = paste0("TMP_lon", my_DE_kriged@sp@coords[,1], "_lat", my_DE_kriged@sp@coords[,2])
tempdat = matrix(test1$val, ncol = length(TMP_var), byrow = T)
kriged_data_wide[,TMP_var] = as.data.frame(tempdat)
# test looks good
image(matrix(unlist(kriged_data_wide[1,2:109]), 
             nrow = length(unique(my_DE_kriged@sp@coords[,1])),
             byrow = F))

# wide variables temperature at lon and lat per day
SLP_var = paste0("SLP_lon", my_DE_kriged_SLP@sp@coords[,1], "_lat", my_DE_kriged_SLP@sp@coords[,2])
slpdat = matrix(test2$val, ncol = length(SLP_var), byrow = T)
kriged_data_wide[,SLP_var] = as.data.frame(slpdat)

saveRDS(kriged_data_wide, "../data/complete_dataset/gridded_wide.RDS")

# ##############################################################################
# create PCA seperately for TMP and SLP
pca_TMP = prcomp(kriged_data_wide[,TMP_var])
plot(pca_TMP)
saveRDS(pca_TMP, "../data/complete_dataset/griddedPC_TMP.RDS")

pc_dat_tmp = data.frame(lon = my_DE_kriged@sp@coords[,1],
                      lat = my_DE_kriged@sp@coords[,2],
                      pca_TMP$rotation)

for (i in 1:10) {
  ggplot() + 
    geom_tile(data = pc_dat_tmp, aes(x=lon, y=lat, fill=eval(parse(text = paste0("PC", i))))) +
    scale_fill_distiller("loading", palette = "Spectral") +
    # geom_point(aes(x = lon, y = lat), data = station_data, size = 1, col = "black") +
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill = NA, colour = "darkgrey", data = ant) +
    coord_cartesian(ylim = c(-90, 0)) +
    theme_bw()
  ggsave(paste0("../plots/covariates/krige_PC_fineGrid/TMP_PC", i, ".png"), height = 5, width = 10)
  # ggsave(paste0("../plots/covariates/krige_PC/TMP_PC", i, ".png"), height = 5, width = 10)
}

# for sea level pressure
pca_SLP = prcomp(kriged_data_wide[,SLP_var])
plot(pca_SLP)
saveRDS(pca_SLP, "../data/complete_dataset/griddedPC_SLP.RDS")

pc_dat_slp = data.frame(lon = my_DE_kriged_SLP@sp@coords[,1],
                        lat = my_DE_kriged_SLP@sp@coords[,2],
                        pca_SLP$rotation)

for (i in 1:10) {
  ggplot() + 
    geom_tile(data = pc_dat_slp, aes(x=lon, y=lat, fill=eval(parse(text = paste0("PC", i))))) +
    scale_fill_distiller("loading", palette = "Spectral") +
    # geom_point(aes(x = lon, y = lat), data = station_data, size = 1, col = "black") +
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill = NA, colour = "darkgrey", data = ant) +
    coord_cartesian(ylim = c(-90, 0)) +
    theme_bw()
  ggsave(paste0("../plots/covariates/krige_PC_fineGrid/SLP_PC", i, ".png"), height = 5, width = 10)
  # ggsave(paste0("../plots/covariates/krige_PC/SLP_PC", i, ".png"), height = 5, width = 10)
}


# ##############################################################################
# create PCA for both TMP and SLP
pca_both = prcomp(kriged_data_wide[,c(TMP_var, SLP_var)])
plot(pca_both)
saveRDS(pca_both, "../data/complete_dataset/griddedPC.RDS")

pc_dat_both = data.frame(lon = my_DE_kriged_SLP@sp@coords[,1],
                        lat = my_DE_kriged_SLP@sp@coords[,2],
                        var = rep(c("TMP", "SLP"), each = length(my_DE_kriged_SLP@sp@coords[,1])),
                        pca_both$rotation)

for (i in 1:5) {
  ggplot() + 
    facet_wrap(~var) +
    geom_tile(data = pc_dat_both, aes(x=lon, y=lat, fill=eval(parse(text = paste0("PC", i))))) +
    scale_fill_distiller("loading", palette = "Spectral") +
    # geom_point(aes(x = lon, y = lat), data = station_data, size = 1, col = "black") +
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill = NA, colour = "darkgrey", data = ant) +
    coord_cartesian(ylim = c(-90, 0)) +
    theme_bw()
  ggsave(paste0("../plots/covariates/krige_PC_fineGrid/Combined_PC", i, ".png"), height = 5, width = 15)
  # ggsave(paste0("../plots/covariates/krige_PC/Combined_PC", i, ".png"), height = 5, width = 15)
}

# ##############################################################################
# # how much efficiency do we gain from grid?
# # compute variance explained by unweighted PCA on grid
# # project gridded data down on original PCA
# # how much efficiency will be gained?

# project gridded data down on original PCA
# get original PCA
cov = readRDS("../data/fogt_predictor_data_imputed/all_monthly_cov.RDS")
# do an eigenvector decomposition on covariates
eigen_cov = prcomp(cov[,-(1:3)])
# variance per PC
plot(eigen_cov$sdev^2)
diag(cov(eigen_cov$x)) %>% plot()
# proportion of variance explained by PC
plot(eigen_cov$sdev^2 / sum(eigen_cov$sdev^2))
(diag(cov(eigen_cov$x)) / sum(diag(cov(eigen_cov$x)))) %>% plot()
(diag(cov(eigen_cov$x)) / sum(diag(cov(eigen_cov$x))))[1:10]
var(eigen_cov$x[,1]) / sum(apply(eigen_cov$x, 2, var))
# divide variance of data from first PC by total sum of variances in data
var(eigen_cov$x[,1]) / sum(apply(cov[,-(1:3)], 2, var))

# ##############################################################################
# predict gridded data from first PC
m_limo = lm(as.matrix(kriged_data_wide[,c(TMP_var, SLP_var)]) ~ eigen_cov$x[,1])
summary(m_limo)

# total variance
sum(diag(cov(kriged_data_wide[,c(TMP_var, SLP_var)])))
# variance explained
sum(diag(cov(fitted(m_limo))))
# residual explained
sum(diag(cov(residuals(m_limo))))
# proportion of variance explained by first PC
sum(diag(cov(fitted(m_limo)))) / sum(diag(cov(kriged_data_wide[,c(TMP_var, SLP_var)])))
# first PC of data explains 9% of variance in gridded data

# ##############################################################################
# # just a sanity check, do it on original data
# predict gridded data from first PC
m_limo_check = lm(as.matrix(cov[,-(1:3)]) ~ eigen_cov$x[,1])
summary(m_limo_check)

# total variance
sum(diag(cov(cov[,-(1:3)])))
# variance explained
sum(diag(cov(fitted(m_limo_check))))
# residual explained
sum(diag(cov(residuals(m_limo_check))))
# proportion of variance explained by first PC
sum(diag(cov(fitted(m_limo_check)))) / sum(diag(cov(cov[,-(1:3)])))
# first PC of data explains 12% of variance in gridded data
# ##############################################################################
# # compute variance explained by original PC
# run loop over first 30 principal components
var_explained = rep(0, 30)
for (i in 1:30) {
  # predict gridded data from second PC
  m_limo = lm(as.matrix(kriged_data_wide[,c(TMP_var, SLP_var)]) ~ eigen_cov$x[,i])
  # proportion of variance explained by first PC
  var_explained[i] = sum(diag(cov(fitted(m_limo)))) / sum(diag(cov(kriged_data_wide[,c(TMP_var, SLP_var)])))
}
plot(var_explained, type = "b")


# #############################################################################
# # variance explained by gridded PC
plot((pca_both$sdev^2 / sum(pca_both$sdev^2))[1:30])
points(var_explained, type = "b", col = 2)

# cumulative variance explained by gridded PC and original PC
# png("../plots/covariates/krige_PC/cumulative_var_explained.png", height = 360, width = 540)
# png("../plots/covariates/krige_PC_fineGrid/cumulative_var_explained.png", height = 360, width = 540)
plot(cumsum((pca_both$sdev^2 / sum(pca_both$sdev^2))[1:30]), 
     type = "b",
     ylim = c(0, 1), 
     xlab = "principal component",
     ylab = "cumulative variance explained")
points(cumsum(var_explained), type = "b", col = 2)
abline(h = c(1), lty = 2)
# dev.off()

# difference between gridded PC and original PC
cumsum((pca_both$sdev^2 / sum(pca_both$sdev^2))[1:30]) - cumsum(var_explained) 



