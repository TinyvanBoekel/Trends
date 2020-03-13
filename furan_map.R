## R file to do regression using map on the furan data at 70 C.
library(ggplot2)
library(rethinking)
library(HDInterval)

time <- c(0,30,60,90,120,150) # in min
furan<- c(0,54,84,147,177,225) # in microgram/L at 70 C
df_furan <- data.frame(time, furan)

# plotting the three priors
prior_c0 <- ggplot(data=data.frame(x=seq(from = -6, to = 10, by = 0.1)), aes(x=x, y=dnorm(x, mean=3, sd=50))) +
  geom_line() +
  ylab("density") + xlab(expression(c[0])) +
  theme_bw()
prior_kr <- ggplot(data=data.frame(x=seq(from = -3, to = 6, by = 0.1)), aes(x=x, y=dnorm(x, mean=1.5, sd=25))) +
  geom_line() +
  ylab("density") + xlab(expression(k[r])) +
  theme_bw()
prior_sigma <- ggplot(data=data.frame(x=seq(from = 0, to = 10, by = 1)), aes(x=x, y=dcauchy(x,location = 3, scale = 25))) +
  geom_line() +
  ylab("density") + xlab(expression(sigma)) +
  theme_bw()

# Multiple plot function, taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

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
multiplot(prior_c0, prior_kr, prior_sigma, cols =3)

#prior predictive simulation
set.seed(200)
N <- 100
c0 <- rlnorm(N,0,1)
kr <- rlnorm(N,1.5,1)
priorsim <- data.frame(c0,kr)


somelines <- ggplot() + 
  geom_abline(data = priorsim, aes(intercept=c0, slope=kr), alpha = 0.5) +
  labs(x = "time (min)", y = "furan") +
  theme_bw()
somelines

# regression using the package rethinking from McElreath

model1_furan <- map(
  alist(
    furan~dnorm(mu,sigma),
    mu <- c0+kr*time,
    c0 ~ dnorm(0,50),
    kr ~ dnorm(1.5,25),
    sigma ~ dcauchy(0,25)
  ),
  data=df_furan, start=list(c0=3, kr=1.5, sigma = 5))

# results of the regression
  precis(model1_furan, corr=TRUE, digits=3) 
  ggplot(data=df_furan, aes(x = time, y = furan))+
    geom_abline(intercept = coef(model1_furan)[1], 
                slope = coef(model1_furan)[2]) +
    geom_point(shape = 1, size = 2, color = "royalblue") +
    theme_bw()

#extracting the posterior  
    post<-extract.samples(model1_furan)
    post
    
    
  # pairs plot with densities
  pairs(model1_furan)
  
# comparison of prior and posterior for each parameter  
  ggplot(post, aes(post$c0))+
    geom_density()+
    geom_line(data=data.frame(x=seq(from = -20, to = 20, by = 0.1)), aes(x=x, y=dnorm(x, mean=0, sd=25)), linetype = 2)+
    geom_vline(xintercept = 3.061, linetype = 3) +
    theme_bw()
  
  ggplot(post, aes(post$kr))+
    geom_density()+
    geom_line(data=data.frame(x=seq(from = -3, to = 6, by = 0.1)), aes(x=x, y=dnorm(x, mean=1.5, sd=25)), linetype = 2)+
    geom_vline(xintercept = 1.485, linetype = 3) +
    theme_bw()
  
  ggplot(post, aes(post$sigma))+
    geom_density()+
    geom_line(data=data.frame(x=seq(from = 0, to = 10, by = 1)), aes(x=x, y=dcauchy(x,location = 3, scale = 20)), linetype = 2)+
    geom_vline(xintercept = 6.328, linetype = 3) +
    theme_bw()
  
# plotting confidence intervals
# density for one particular time value (50 in this case)
  mu_50 <- data.frame(post$c0 + post$kr*50)
  ggplot(mu_50) +
    geom_density(aes(x=post$c0+post$kr*50))+
    theme_bw() + 
    labs(x=expression(mu[50]), caption = expression(paste("Density of ", mu[50])))
  
#plotting hdpi interval:

  hdi(mu_50, credMass = 0.95)
  ggplot(mu_50) +
    geom_density(aes(x=post$c0+post$kr*50))+
    geom_vline(xintercept = hdi(mu_50, credMass = .95),
               color = "royalblue4", linetype = 2) +
    theme_bw() + 
    labs(x=expression(mu[50]), caption = expression(paste("Density of ", mu[50])))
  
# plotting values of the mean over the whole time period
  time.seq <- seq(from =0, to =150, by=1)
  mu<-link(model1_furan, data=data.frame(time=time.seq))
  mu.mean <- apply(mu, 2, mean)
  mu.HPDI <-  apply(mu, 2, HPDI, prob = 0.95)
#  plot(furan~time, xlab=("time (min)"), ylab=("furan (mg/L)"),df_furan)
#  lines(time.seq, mu.mean)
#  shade(mu.HPDI, time.seq)
 
 # plotting prediction intervals
  sim.furan <- sim(model1_furan, data = list(time=time.seq))
  furan.PI <- apply(sim.furan, 2, PI, prob = 0.95)
  plot(furan~time, xlab=("time (min)"), ylab=("furan (microgram/L)"),df_furan)
  lines(time.seq, mu.mean)
  shade(mu.HPDI, time.seq)
  shade(furan.PI, time.seq)
  
  