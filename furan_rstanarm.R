# script on how to use Rstanarm for the furan case: zero order reaction

library(rstanarm)
library("ggplot2")
library("broom")
library("bayesplot")
library(dplyr)
library(coda)

#data
time <- c(0,30,60,90,120,150) # in min
furan<- c(0,54,84,147,177,225) # in microgram/L at 70 C
df_furan <- data.frame(time, furan)

furan_fit <- stan_glm(furan~time, data=df_furan, iter=4000, warmup=2000, chains =4, 
                           prior_intercept = normal(0,50), prior = normal(0,25),
                           prior_aux = cauchy(0,25))

plot(furan_fit, "trace") #check the trace plots

summary(furan_fit) #all possible information together

coef(furan_fit) #shows only the parameters

# select samples from the posterior to draw the regression line

post <- furan_fit %>%   as_tibble %>% rename(intercept = `(Intercept)`) %>% 
  select(-sigma)

post # tibble that contains 8000 samples of intercepts and slopes

#plot 1000 out of these 8000 samples
n_draws <- 1000
alpha_level <- .15
col_draw <- "grey60"
col_median <-  "#3366FF"

ggplot(df_furan) + 
  aes(x = time, y = furan) + 
  geom_abline(aes(intercept = intercept, slope = time), 
              data = sample_n(post, n_draws), color = col_draw, 
              alpha = alpha_level) + 
  geom_abline(intercept = median(post$intercept), 
              slope = median(post$time), 
              size = 1, color = col_median) +
  geom_point()


posterior <- as.array(furan_fit)
#dim(posterior)
#dimnames(posterior)
color_scheme_set("red")
mcmc_intervals(posterior, pars = c("(Intercept)", "time", "sigma"))
mcmc_areas(
  posterior,
  pars = c("(Intercept)", "time", "sigma"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

color_scheme_set("green")
mcmc_hist(posterior, pars = c("(Intercept)", "time","sigma"))

color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("(Intercept)", "time", "sigma"),
           off_diag_args = list(size = 1.5))


mcmc_combo(furan_fit)


#regression line with 95% prediction interval

newdata = data.frame(time = seq(min(df_furan$time, na.rm = TRUE), max(df_furan$time, 
                                                                na.rm = TRUE),len = 1000))
fit = posterior_predict(furan_fit, newdata = newdata)

newdata = newdata %>% cbind(tidyMCMC(as.mcmc(fit), conf.int = TRUE, conf.method = "HPDinterval"))

ggplot(newdata, aes(y = estimate, x = time)) + geom_point(data = df_furan, aes(y = furan))+ 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "blue", alpha = 0.3) + 
  scale_y_continuous("furan") +
  scale_x_continuous("time") + 
  theme_classic()

pp_check(furan_fit)

# a very nice way to analyse the regression results is via shinystan:

launch_shinystan(furan_fit)
