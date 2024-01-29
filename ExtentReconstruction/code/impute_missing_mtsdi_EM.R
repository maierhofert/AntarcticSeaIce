library("mgcv")
library("mtsdi")

# detrend by smooth day of year effect
gamTMAX = mgcv::gam(TMAX ~ s(Doy, bs = "cc"), data = datASN2011)
datASN2011$resTMAX = datASN2011$TMAX - predict(gamTMAX, newdata = datASN2011)
gamTMIN = mgcv::gam(TMIN ~ s(Doy, bs = "cc"), data = datASN2011)
datASN2011$resTMIN = datASN2011$TMIN - predict(gamTMIN, newdata = datASN2011)

# scale the residuals
datASN2011$sTMAX = scale(datASN2011$resTMAX, center = TRUE, scale = TRUE)
datASN2011$sTMIN = scale(datASN2011$resTMIN, center = TRUE, scale = TRUE)

# add linear time
datASN2011$ts = 1:nrow(datASN2011)

# # fit imputation model
# tmod = mnimput(formula = ~ sTMIN + sTMAX, 
#                # dataset = datASN2011,
#                dataset = datASN2011[450:550,],
#                # dataset = datASN2011[1:1000,], 
#                # dataset = datASN2011[22000:23000,], 
#                ts = TRUE,
#                method="spline", 
#                sp.control=list(df=5))

# fit imputation model
# with annual cycle included in state space gam
tmod = mnimput(formula = ~ TMIN + TMAX, 
               dataset = datASN2011,
               # dataset = datASN2011[450:550,],
               # dataset = datASN2011[1:1000,], 
               # dataset = datASN2011[22000:23000,], 
               ts = TRUE,
               method="gam", 
               ga.control=list(formula = c(TMIN ~ s(ts) + s(Doy), 
                                           TMAX ~ s(ts) + s(Doy))))
summary(tmod)
# plot(tmod)
plot(tmod, xlim = c(450, 550))
plot(tmod, xlim = c(1, 1000))

# debugonce(mtsdi::em.recursion)

# # hack the em.gam function to use mgcv::gam
# my_em.gam = function (formula, xn, names, dataset, w, eps, maxit, bf.eps, bf.maxit) {
#   rows <- dim(xn)[1]
#   cols <- dim(xn)[2]
#   t <- seq(1:rows)
#   ga.pred <- matrix(NA, nrow = rows, ncol = cols)
#   if (is.null(w)) 
#     w <- rep(1, rows)
#   models <- as.list(rep(NA, cols))
#   XN <- as.data.frame(xn)
#   dimnames(XN) <- names
#   dataset <- cbind.data.frame(xn, dataset)
#   dataset$w = w
#   # formula <- paste("XN$", formula, sep = "")
#   for (j in 1:cols) {
#     # models[[j]] <- mgcv::gam(as.formula(formula[[j]]), family = gaussian(), 
#     #                          data = dataset, weights = w, na.action = na.exclude)
#     models[[j]] <- gam(as.formula(formula[[j]]), family = gaussian(), 
#                              data = dataset, weights = w, na.action = na.exclude)
#     
#     # epsilon = eps, maxit = maxit, bf.epsilon = bf.eps, 
#     # bf.maxit = bf.maxit)
#     ga.pred[, j] <- fitted(models[[j]])
#     # ga.pred[, j] <- fitted(models[[j]])
#   }
#   retval <- list(ga.pred = ga.pred, models = models)
#   return(retval)
# }
# assignInNamespace("em.gam", my_em.gam, "mtsdi")



# fit imputation model
# with detrended values
tmod = mnimput(formula = ~ sTMIN + sTMAX, 
               # dataset = datASN2011,
               # dataset = datASN2011[450:550,],
               # dataset = datASN2011[1:1000,], 
               dataset = datASN2011[22000:23000,], 
               ts = TRUE,
               method="gam", 
               ga.control=list(formula = c(sTMIN~s(ts, df = 10), 
                                           sTMAX~s(ts, df = 10))))
summary(tmod)
plot(tmod, xlim = c(450, 550))
plot(tmod, xlim = c(1, 1000))
# imputed values seem overly smooth, see TMIN around 900

# add imputed values to data
datASN2011[, c("predTMIN", "predTMAX")] = predict(tmod)


