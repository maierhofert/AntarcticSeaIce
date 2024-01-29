library("dplyr")
library(rstan)
rstan_options(auto_write = TRUE)


library("brms")
# devtools::install_github("paul-buerkner/brms", dependencies = TRUE)

options(mc.cores = 10)

# get data
stan_data = list(K = 5,
                 J = 10,
                 N = 501,
                 x = all_dat[11000:11500,200:209],
                 y = all_dat[11000:11500,4:8])


# set up multivariate arma model
# is ar good enough?
this.formula = brmsformula(
  mvbind(King_Hakon, Ross, East_Antarctica, Weddell, Bellingshausen_Amundsen_Sea) ~ VMarsh_O_Higgins_tmp + VScott_Base_mcMurdo_SLP + VScott_Base_McMurdo_tmp,
  autocor = ~ arma(p = 1, q = 1))
fit1 <- brm(
  this.formula,
  data = all_dat[11000:11500,], chains = 2, cores = 2, iter = 200
)
# change adapt parameter
summary(fit1)


# set up multivariate ar model with horseshoe prior
this.formula2 = brmsformula(
  mvbind(King_Hakon, Ross, East_Antarctica, Weddell, Bellingshausen_Amundsen_Sea) ~ VMarsh_O_Higgins_tmp + VScott_Base_mcMurdo_SLP + VScott_Base_McMurdo_tmp,
  autocor = ~ ar(p = 1))
fit2 <- brm(
  this.formula2,
  prior = set_prior(horseshoe(df = 1, par_ratio = 0.1)),
  data = all_dat[11000:11500,], chains = 2, cores = 2, iter = 200
)
summary(fit2)

# fit2 <- brm(
#   this.formula2,
#   prior = hs(df = 1, global_scale = 0.1),
#   data = all_dat[11000:11500,], chains = 2, cores = 2, iter = 200
# )
# summary(fit2)



# set up multivariate ar model with seasonally varying effects
this.formula3 = brmsformula(
  mvbind(King_Hakon, Ross) ~ -1 + s(Doy,by=VMarsh_O_Higgins_tmp,k=4,bs="cc"),
  autocor = ~ ar(p = 1))

this.formula3 = brmsformula(
  mvbind(King_Hakon, Ross) ~ -1 + s(Doy, by = VMarsh_O_Higgins_tmp, k = 4, bs = 'cc'),
  autocor = ~ ar(p = 1))

this.prior = get_prior(this.formula3, data = all_dat[11000:11500,])
# reparametrize 1 / (1 + var) limited to [0, 1]?
this.prior = c(# prior(normal(100, 10), class = Intercept),
  # prior(student_t(3, 0, 1), class = ar, resp = KingHakon),
  prior(uniform(0, 1), class = ar, resp = KingHakon),
  prior(cauchy(0, 0.001), class = b, resp = KingHakon),
  prior(cauchy(0, 0.001), class = sds, resp = KingHakon),
  prior(exponential(1), class = sigma, resp = KingHakon),
  #
  prior(uniform(0, 1), class = ar, resp = Ross),
  prior(cauchy(0, 0.001), class = b, resp = Ross),
  prior(cauchy(0, 0.001), class = sds, resp = Ross),
  prior(exponential(1), class = sigma, resp = Ross))

# this.prior = get_prior(this.formula3, data = all_dat[11000:11500,])
# this.prior
# str(this.prior)
# this.prior$prior[5] = "student_t(3, 0, 0.5)"
# this.prior$prior[6] = "cauchy(0, 0.1)"

this.prior

fit3 <- brm(
  this.formula3,
  prior = this.prior,
  data = all_dat[11000:11500,], chains = 2, cores = 2, iter = 200
)
# brms::prior_summary(fit3)
summary(fit3)
# stancode(fit3)

plot(fit3)
str(fit3)

# marg_smooth = marginal_smooths(fit3)
# plot(marg_smooth)

# debugonce(brms:::conditional_smooths.mvbrmsterms)
# debugonce(brms:::conditional_smooths.brmsterms)
# debugonce(brms:::conditional_smooths.btl)
# debugonce(brms:::posterior_smooths.btl)

# you need to hack into the brms package here
# source("my_posterior_smooths.R")
cond_smooth = conditional_smooths(fit3)
plot(cond_smooth)

library("ggplot2")

# plot estimates
predict(fit3, newdata = all_dat[12500:13000,-c(4:8)])[,,"KingHakon"] %>% 
  data.frame() %>% 
  bind_cols(all_dat[12500:13000,]) %>% 
  ggplot(aes(x = Date, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  # geom_hline(yintercept = fixef(b4.11)[1, 1], color = "white", linetype = 2) +
  geom_ribbon(fill = "grey83") +
  geom_line() +
  geom_line(aes(x = Date, y = King_Hakon), color = "blue") +
  theme_minimal()

predict(fit3, newdata = all_dat[12500:13000,-c(4:8)])[,,"Ross"] %>% 
  data.frame() %>% 
  bind_cols(all_dat[12500:13000,]) %>% 
  ggplot(aes(x = Date, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  # geom_hline(yintercept = fixef(b4.11)[1, 1], color = "white", linetype = 2) +
  geom_ribbon(fill = "grey83") +
  geom_line() +
  geom_line(aes(x = Date, y = Ross), color = "blue") +
  theme_minimal()

fixef(fit3)



