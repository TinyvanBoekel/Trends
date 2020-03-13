## R file to do regression using map on the furan data at 70 C.
library(ggplot2)
library(brms)
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


#regression using brms
furanfit <- 
  brm(data = furan_data, family = gaussian,
      formula = furan ~ 1 + time,
      prior = c(set_prior("normal(0, 50)", class = "Intercept"),
                set_prior("normal(0, 25)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, 
      control = list(adapt_delta =0.95), file="furanfit")

# results of the regression
plot(furanfit, type="trace") 
summary(furanfit)
 
#extracting the posterior  
furanpost <- posterior_samples(furanfit)
    
  # pairs plot with densities
pairs(furanfit)
  
