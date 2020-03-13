library(papaja)
library(rmarkdown)
library(bookdown)
library(ggplot2)
library(dplyr)
library(brms)
library(tidyverse)
library(rstan)
library(broom)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

#algorithm needed to plot figures next to each other
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#preparing Figure 2
#normal

plot_normal <- ggplot(data.frame((x=c(0,100))), aes(x=x))+
  stat_function(fun=dnorm, args=list(50,10), geom = "area", alpha = 1,
                fill = "lightblue")+
  stat_function(
    fun = dnorm, args = list(50,10),
    geom = "line", alpha = 1)+
  theme_bw()+
  labs(x="x", y="density")+
  annotate("text", x=100, y=0.04, label="A")

# uniform
plot_uniform <- ggplot(data.frame((x=c(0,100))), aes(x=x))+
  stat_function(fun=dunif, args=list(10,50), geom="area", alpha=1,
                fill="lightblue")+
  stat_function(fun=dunif, args=list(10,50), geom="line", alpha=1)+
  theme_bw()+
  labs(x="x", y="density")+
  annotate("text", x=100, y=0.025, label="B")

# gamma

plot_gamma <- ggplot(data.frame((x=c(0,100))), aes(x=x))+
  stat_function(fun=dgamma, args=list(2,0.1), geom="area", alpha=1,
                fill="lightblue")+
  stat_function(fun=dgamma, args=list(2,0.1), geom="line", alpha=1)+
  theme_bw()+
  labs(x="x", y="density")+
  annotate("text", x=100, y=0.035, label="C")

# half-cauchy
plot_hcauchy <- ggplot(data.frame((x=c(0,100))), aes(x=x))+
  stat_function(fun=dcauchy, args=list(0,25), geom="area", alpha=1,
                fill="lightblue")+
  stat_function(fun=dcauchy, args=list(0,25), geom="line", alpha=1)+
  theme_bw()+
  labs(x="x", y="density")+
  annotate("text", x=100, y=0.0125, label="D")

#Figure 2
multiplot(plot_normal, plot_gamma, plot_uniform, plot_hcauchy, cols=2)

ndata <- 20
theta.data <- 12
variance.known <- 7
prior.theta <- 5
prior.variance1 <- 50
prior.sd1 <- sqrt(prior.variance1)

x <- rnorm(ndata, mean=theta.data, sd=sqrt(variance.known))
xlim <- theta.data+c(-1, 1)*4*sqrt(variance.known)
steplength <- (xlim[2]-xlim[1])/500
xx <- seq(xlim[1], xlim[2], by=steplength)

#calculate the likelihood function

likelihood.function <- function(x, theta, variance.known) prod(1/(sqrt(2*pi*variance.known))*exp(-1/(2*variance.known)*(x-theta)^2))

likelihood <- numeric(length(xx))
for(i in 1:length(xx)) likelihood[i] <- likelihood.function(x, xx[i], variance.known=variance.known)

lf <- likelihood.function
lf2 <- likelihood/sum(likelihood)*1/steplength

lf_data <- data.frame(xx, lf2)

#calculate the posterior

posterior.mean1 <- (prior.theta/prior.variance1 + ndata*mean(x)/variance.known)/(1/prior.variance1 + ndata/variance.known)
posterior.variance1 <- 1/(1/prior.variance1 + ndata/variance.known)
posterior.sd1 <- sqrt(posterior.variance1)

plot1 <- ggplot(data.frame(x=xx), aes(x=xx))+
  stat_function(fun=dnorm, args=list(mean=posterior.mean1, sd=posterior.sd1))+
  stat_function(fun = dnorm, args=list(mean=prior.theta, sd=prior.sd1), linetype=2)+
  theme_bw()+labs(x=expression(theta), y="density")+
  geom_line(data=lf_data, aes(x=lf_data$xx, y=lf_data$lf2), linetype=3)+
  geom_text(x=18, y=0.64, label=expression(sigma^{2}==50))

#other value of prior variance

prior.variance2 <- 5
prior.sd2 <- sqrt(prior.variance2)

#calculate the posterior

posterior.mean2 <- (prior.theta/prior.variance2 + ndata*mean(x)/variance.known)/(1/prior.variance2 + ndata/variance.known)
posterior.variance2 <- 1/(1/prior.variance2 + ndata/variance.known)
posterior.sd2 <- sqrt(posterior.variance2)

plot2 <- ggplot(data.frame(x=xx), aes(x=xx))+
  stat_function(fun=dnorm, args=list(mean=posterior.mean2, sd=posterior.sd2))+
  stat_function(fun = dnorm, args=list(mean=prior.theta, sd=prior.sd2), linetype=2)+
  theme_bw()+labs(x=expression(theta), y="density")+
  geom_line(data=lf_data, aes(x=lf_data$xx, y=lf_data$lf2), linetype=3)+
  geom_text(x=18, y=0.65, label=expression(sigma^{2}==5))

# other value of prior variance
prior.variance3 <- 0.5
prior.sd3 <- sqrt(prior.variance3)

#calculate the posterior

posterior.mean3 <- (prior.theta/prior.variance3 + ndata*mean(x)/variance.known)/(1/prior.variance3 + ndata/variance.known)
posterior.variance3 <- 1/(1/prior.variance3 + ndata/variance.known)
posterior.sd3 <- sqrt(posterior.variance3)

plot3 <- ggplot(data.frame(x=xx), aes(x=xx))+
  stat_function(fun=dnorm, args=list(mean=posterior.mean3, sd=posterior.sd3))+
  stat_function(fun = dnorm, args=list(mean=prior.theta, sd=prior.sd3), linetype=2)+
  theme_bw()+labs(x=expression(theta), y="density")+
  geom_line(data=lf_data, aes(x=lf_data$xx, y=lf_data$lf2), linetype=3)+
  geom_text(x=17, y=0.82, label=expression(sigma^{2}==0.5))

ndata4 <- 20
theta.data <- 12
variance.known <- 7
prior.theta <- 5
prior.variance <- 0.5
prior.sd <- sqrt(prior.variance)

x <- rnorm(ndata4, mean=theta.data, sd=sqrt(variance.known))
xlim <- theta.data+c(-1, 1)*4*sqrt(variance.known)
steplength <- (xlim[2]-xlim[1])/500
xx <- seq(xlim[1], xlim[2], by=steplength)

#calculate the likelihood function

likelihood.function <- function(x, theta, variance.known) prod(1/(sqrt(2*pi*variance.known))*exp(-1/(2*variance.known)*(x-theta)^2))

likelihood <- numeric(length(xx))
for(i in 1:length(xx)) likelihood[i] <- likelihood.function(x, xx[i], variance.known=variance.known)

lf <- likelihood.function
lf2 <- likelihood/sum(likelihood)*1/steplength

lf_data <- data.frame(xx, lf2)

#calculate the posterior

posterior.mean4 <- (prior.theta/prior.variance + ndata4*mean(x)/variance.known)/(1/prior.variance + ndata4/variance.known)
posterior.variance4 <- 1/(1/prior.variance + ndata4/variance.known)
posterior.sd4 <- sqrt(posterior.variance4)

plot4 <- ggplot(data.frame(x=xx), aes(x=xx))+
  stat_function(fun=dnorm, args=list(mean=posterior.mean4, sd=posterior.sd4))+
  stat_function(fun = dnorm, args=list(mean=prior.theta, sd=prior.sd), linetype=2)+
  theme_bw()+labs(x=expression(theta), y="density")+
  geom_line(data=lf_data, aes(x=lf_data$xx, y=lf_data$lf2), linetype=3)+
  geom_text(x=17, y=0.82, label="n=20")

#other value of number of data

ndata5 <- 50

#calculate the posterior

posterior.mean5 <- (prior.theta/prior.variance + ndata5*mean(x)/variance.known)/(1/prior.variance + ndata5/variance.known)
posterior.variance5 <- 1/(1/prior.variance + ndata5/variance.known)
posterior.sd5 <- sqrt(posterior.variance5)

plot5 <- ggplot(data.frame(x=xx), aes(x=xx))+
  stat_function(fun=dnorm, args=list(mean=posterior.mean5, sd=posterior.sd5))+
  stat_function(fun = dnorm, args=list(mean=prior.theta, sd=prior.sd), linetype=2)+
  theme_bw()+labs(x=expression(theta), y="density")+
  geom_line(data=lf_data, aes(x=lf_data$xx, y=lf_data$lf2), linetype=3)+
  geom_text(x=17, y=1.12, label="n=50")

# other value of number of data
ndata6 <- 100

#calculate the posterior

posterior.mean6 <- (prior.theta/prior.variance + ndata6*mean(x)/variance.known)/(1/prior.variance + ndata6/variance.known)
posterior.variance6 <- 1/(1/prior.variance + ndata6/variance.known)
posterior.sd6 <- sqrt(posterior.variance6)

plot6 <- ggplot(data.frame(x=xx), aes(x=xx))+
  stat_function(fun=dnorm, args=list(mean=posterior.mean6, sd=posterior.sd6))+
  stat_function(fun = dnorm, args=list(mean=prior.theta, sd=prior.sd), linetype=2)+
  theme_bw()+labs(x=expression(theta), y="density")+
  geom_line(data=lf_data, aes(x=lf_data$xx, y=lf_data$lf2), linetype=3)+
  geom_text(x=18, y=1.4, label="n=100")

#plotting Figure 3
multiplot(plot1,plot4,plot2, plot5, plot3, plot6, cols = 3)


#furan data, source: Huang & Barringer, LWT/Food Sci Technol. 67(2016)200-205

time <- c(0,30,60,90,120,150) # in min
furan<- c(0,54,84,147,177,225) # in microgram/L
furandata <- data.frame(time, furan)

# doing the Bayesian regression with brms
furanfit <- 
  brm(data = furandata, family = gaussian,
      formula = furan ~ 1 + time,
      prior = c(set_prior("normal(0, 100)", class = "Intercept"),
                set_prior("normal(1, 10)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      control = list(adapt_delta =0.95), 
      inits = 0, seed=15, file="furanfit")

#extracting the posterior
furanpost <- posterior_samples(furanfit)

#preparing Figure 4
lines10 <- furanpost[1:10, 1:2]
lines10$iternum <- seq(1, length(lines10$b_Intercept))
somelines <- ggplot(data= furandata, aes(x=time, y = furan)) + 
  geom_point(size =2) + 
  geom_abline(data = lines10, aes(intercept=b_Intercept, slope=b_time, group = iternum), alpha = 0.5) +
  labs(x = "time (min)", y = "furan") +
  theme_bw() +
 annotate("text", x=10, y=200, label="A")

# calculating the regression line: 
regrline <-  ggplot(data=furandata, aes(x=time, y=furan))+
  geom_point(size = 2) +
  geom_abline(data=furanpost, aes(intercept=mean(b_Intercept), slope=mean(b_time)))+
  theme_bw() +
  labs(x = "time (min)", y = "furan") +
annotate("text", x=10, y=200, label="B")

#Figure 4
multiplot(somelines,regrline, cols=2) 

#preparing Figure 5

mu_100 <- data.frame(furanpost$b_Intercept + furanpost$b_time*100)
mu2_100 <- rnorm(n=nrow(furanpost), mean=furanpost$b_Intercept+furanpost$b_time*100, sd=furanpost$sigma)
plotA <- ggplot(mu_100) +
  geom_density(aes(x=furanpost$b_Intercept+furanpost$b_time*100))+
  geom_density(aes(x=rnorm(n=nrow(furanpost), mean=furanpost$b_Intercept+furanpost$b_time*100, sd=furanpost$sigma)), lty=2)+
  theme_bw() + 
  labs(x=expression(mu[100]))+
  coord_cartesian(xlim=c(120,180))+ 
  geom_vline(xintercept = mean(furanpost$b_Intercept)+mean(furanpost$b_time)*100, lty = 4)+
  annotate("text", x=125, y=0.09, label="A")

time.seq <- data.frame(time = seq(from = 0, to = 160, by = 1))

#fitted is about mu:
muSummary <-
  fitted(furanfit, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

#predict is about future individual values:
pred.furan <-
  predict(furanfit,
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

#plot of fitted and predicted values
plotB <- furandata %>%
  ggplot(aes(x = time, y = furan)) +
  geom_ribbon(data = pred.furan, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "grey83") +
  geom_ribbon(data = muSummary, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "grey70") +
  geom_line(data = muSummary, aes(y = Estimate)) +
  geom_point(color = "navyblue", shape = 16, size = 1.5, alpha = 2/3) +
  theme_bw()+
  annotate('text', x=5, y=260, label="B")

#Figure 5
multiplot(plotA, plotB, cols = 2)

#plotting priors and posteriors

c0 <- ggplot(data=furanpost, aes(x=b_Intercept))+
  geom_density()+
  geom_line(data=data.frame(x=seq(from = -100, to = 100, by = 0.1)), aes(x=x, y=dnorm(x, mean=0, sd=100)), linetype = 2)+
  tidybayes::stat_pointintervalh(aes(y=0), .width = 0.95)+
  theme_bw()+
  labs(x = expression(c[0]), y = "density")+
  theme_bw()+ 
  annotate("text", x=75, y=0.06, label="A")

kr <- ggplot(data=furanpost, aes(x=b_time))+
  geom_density()+
  geom_line(data=data.frame(x=seq(from = 0, to = 3, by = 0.01)), aes(x=x, y=dnorm(x, mean=1, sd=10)), linetype = 2)+
  tidybayes::stat_pointintervalh(aes(y=0), .width = 0.95)+
  theme_bw()+
  labs(x = expression(k[r]), y = "density")+
  theme_bw()+
  annotate("text", x=2.5, y=5.5, label="B")

sigma <- ggplot(data=furanpost, aes(x=sigma))+
  geom_density()+
  geom_line(data=data.frame(x=seq(from = 0, to = 75, by = 0.1)), aes(x=x, y=dcauchy(x,location = 0, scale = 25)), linetype = 2)+
  theme_bw() +
  labs(x = expression(sigma), y = "density")+
  tidybayes::stat_pointintervalh(aes(y=0), .width = 0.95)+
  theme_bw()+
  annotate("text", x=70, y=0.12, label="C")

#Plotting Figure 6
multiplot(c0, kr, sigma, cols=3)

# Figures and Tables in the supplement
#Table 1 

(ref:Tab1) Numerical output of the function 'lm' for linear regression of the zero-order model for the formation of furan in soy sauce at 70 ^0^C.

```{r linregr, echo=FALSE, warning=FALSE, results="asis"}
result_lm <- lm(furan~1+time, data=furandata)

apa_lm <- apa_print(result_lm)

apa_lm$table %>% rename(parameter=predictor)
names(apa_lm$table)[names(apa_lm$table)=="b"] <- "estimate"

#write.csv(summary_mod1, file = "Table_2.csv", row.names = FALSE)

apa_table(apa_lm$table, caption="(ref:Tab1)")
```

#Table 2 here.

(ref:Tab2) Numerical results for the Bayesian estimation of parameters in the zero-order model describing formation of furan in soy sauce at 70 ^0^C.

```{r numsummary, echo=FALSE, warning=FALSE, results="asis"}

a <- summary(furanfit)
summary_mod1 <- rbind(data.frame(a$fixed), data.frame(a$spec_pars) )
rownames(summary_mod1) <- c("$c_0$", "$k_r$", "$\\sigma$")
colnames(summary_mod1) <- c("mean","SE", "lower bound", "upper bound", "ESS", "Rhat")

#summary_mod1 %>% rownames_to_column(var = "parameter")

write.csv(summary_mod1, file = "Table_2a.csv", row.names = FALSE)

apa_table(
  summary_mod1,
  placement = "H",
  align = c("c", "c", "c", "c", "c", "c"),
  caption = "(ref:Tab2)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,2,2,2,0,3)
  ),
  escape = FALSE
)
```
#Table 3.  
(ref:Tab3) The first 6 entries from the posterior distribution obtained with brms for the zero-order model, showing estimates of the intercept ($c_0$), the slope ($k_r$), and standard deviation sigma. 

```{r furanpost, echo=FALSE, results="asis"}

sample_data <- head(furanpost)
sample_data <- sample_data %>% select(-lp__) %>% mutate(Entry=c(1,2,3,4,5,6)) %>% select(Entry, b_Intercept,b_time,sigma)

colnames(sample_data) <- c("Entry","intercept", "slope", "sigma")
write.csv(sample_data, file = "Table_3a.csv", row.names = FALSE)

apa_table(
  sample_data,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:Tab3)",
  note = NULL,
  small = TRUE,
  format.args = list(
    digits = c(0,2,2,2),
    decimal.mark = ".", big.mark = ""
  )
)
```
#Table 4
(ref:Tab4) Estimates for the zero-order parameters of the furan regression model using the Student distribution for the likelihood.

```{r numsummary3, echo=FALSE, warning=FALSE, results="asis"}

cc <- summary(furanfit_student)
summary_mod3 <- rbind(data.frame(cc$fixed), data.frame(cc$spec_pars) )
rownames(summary_mod3) <- c("$c_0$", "$k_r$", "$\\sigma$", "$\\nu$")
colnames(summary_mod3) <- c("mean","SE", "lower bound", "upper bound", "ESS", "Rhat")
#summary_mod3 %<>% rownames_to_column(var = "parameter")
write.csv(summary_mod3, file = "Table_4a.csv", row.names = FALSE)

apa_table(
  summary_mod3,
  placement = "H",
  align = c("c", "c", "c", "c", "c", "c"),
  caption = "(ref:Tab4)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(0,2,2,2,2,0,3)
  ),
  escape = FALSE
)
```

#preparing Figure S1

trace1 <- stanplot(furanfit, type="trace")+
  labs(x = "no. of iterations after warmup", y = "") 
library(GGally)
furanpost1 <- posterior_samples(furanfit)

corplot1 <- furanpost1 %>% select(b_Intercept:sigma) %>% ggpairs(diag=list(continuous="barDiag")) + theme_bw()

#FPlotting igure S1
multiplot(trace1, corplot1, cols = 1)

#doing Bayesian regression with Student distribution
furanfit_student <- 
  brm(data = furandata, family = student,
      formula = furan ~ 1 + time,
      prior = c(set_prior("normal(0, 100)", class = "Intercept"),
                set_prior("normal(1, 10)", class = "b"),
                set_prior("gamma(2,0.1)", class = "nu"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      control = list(adapt_delta =0.95),
      seed=15, file="furanfit_student")

#preparing Figure S2
trace2 <- stanplot(furanfit_student, type="trace")+
  labs(x = "no. of iterations after warmup", y = "") 

furanpost2 <- posterior_samples(furanfit_student)

corplot2 <- furanpost2 %>% select(b_Intercept:nu) %>% ggpairs(diag=list(continuous="barDiag")) + theme_bw()

#Plotting Figure S2
multiplot(trace2, corplot2, cols = 1)
