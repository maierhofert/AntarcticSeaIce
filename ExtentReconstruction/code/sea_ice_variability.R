
library("dplyr")
raw_dat = read.csv("../data/nsidc_sea_ice/S_seaice_extent_daily.csv")
raw_dat = raw_dat[-1, -6]
raw_dat = apply(raw_dat, 2, as.numeric) %>% as.data.frame()

raw_dat_mon = raw_dat %>%
  group_by(Year, Month) %>%
  summarize(Extent = mean(Extent))
raw_dat_mon = raw_dat_mon[-(1:3),]
tail(raw_dat_mon, 8)

# compute
monthly_mean = raw_dat_mon %>%
  subset(Year <= 2020) %>%
  group_by(Month) %>%
  summarize(Extent = mean(Extent))

# compute anomaly
raw_dat_mon$Anomaly = raw_dat_mon$Extent - c(rep(monthly_mean$Extent, 43),
                                             monthly_mean$Extent[1:8])

tail(raw_dat_mon, 8)
raw_dat_mon$tdat = raw_dat_mon$Year + (raw_dat_mon$Month - 0.5) / 12

# anomaly and extent
plot(raw_dat_mon$tdat, raw_dat_mon$Anomaly)
plot(raw_dat_mon$tdat, raw_dat_mon$Extent, type = "l")
#
plot(raw_dat_mon$tdat[493:524], raw_dat_mon$Anomaly[493:524])
plot(raw_dat_mon$tdat[493:524], raw_dat_mon$Extent[493:524])

tail(raw_dat_mon, 8)

# the variability varies from month to month
# variability in some months is bout twice as high as in others
raw_dat_mon %>%
  group_by(Month) %>%
  summarize(sd_Extent = sd(Extent)) %>%
  plot(type = "l", ylim = c(0,1.3), col = 2)
# raw_dat_mon %>%
#   +   group_by(Month) %>%
#   +   summarize(sd_Extent = sd(Extent))
# # A tibble: 12 × 2
# Month sd_Extent
# <dbl>     <dbl>
# 1     1     0.699
# 2     2     0.431
# 3     3     0.570
# 4     4     0.712
# 5     5     0.671
# 6     6     0.598
# 7     7     0.461
# 8     8     0.435
# 9     9     0.421
# 10    10    0.401
# 11    11    0.503
# 12    12    0.921

raw_dat_mon %>%
  group_by(Month) %>%
  summarize(Extent = mean(Extent) / 15 ) %>%
  lines()

# on log scale the variaibility of low sea ice months is 7 times higher than
# high sea ice months
raw_dat_mon %>%
  mutate(logExtent = log(Extent)) %>%
  group_by(Month) %>%
  summarize(sd_logExtent = sd(logExtent)) %>%
  plot(type = "l")

# on sqrt scale the variaibility of low sea ice months is 7 times higher than
# high sea ice months
raw_dat_mon %>%
  mutate(sqrtExtent = sqrt(Extent)) %>%
  group_by(Month) %>%
  summarize(sd_sqrtExtent = sd(sqrtExtent)) %>%
  plot(type = "l")


# nthroot = function(x, n){x^(1/n)}
# raw_dat_mon %>%
#   mutate(nthrootExtent = nthroot(Extent, n = 2)) %>%
#   group_by(Month) %>%
#   summarize(sd_nthrootExtent = sd(nthrootExtent)) %>%
#   plot(type = "l")
#
# exp_b = function(x, b){b^x}
# raw_dat_mon %>%
#   mutate(exp_bExtent = exp_b(Extent, b = 1.1)) %>%
#   group_by(Month) %>%
#   summarize(sd_exp_bExtent = sd(exp_bExtent)) %>%
#   plot(type = "l")
#
# raw_dat_mon %>%
#   mutate(logbExtent = log(Extent, base = 3)) %>%
#   group_by(Month) %>%
#   summarize(sd_logbExtent = sd(logbExtent)) %>%
#   plot(type = "l")


# plot monthly sd and monthly squared derivative, these line up very nicely
raw_dat_mon %>%
  group_by(Month) %>%
  summarize(sd_Extent = sd(Extent)) %>%
  plot(type = "l", ylim = c(0,1), col = 2)

monthly_mean = raw_dat_mon %>%
  group_by(Month) %>%
  summarize(Extent = mean(Extent) / 15 )
monthly_der = c(diff((monthly_mean$Extent)),
                monthly_mean$Extent[12] - monthly_mean$Extent[1])
monthly_der_sq = monthly_der^2
lines(monthly_der_sq*5)
