# linear regression model for furan data with Stan
# after the tutorial stan-intro from ourcodingclub
library(rstan)

# data:
time <- c(0,30,60,90,120,150) # in min
furan<- c(0,54,84,147,177,225) # in microgram/L
N<-length(time)
stan_data <- list(N=N, x=time, y = furan)

# stancode, format write("code","filename"):

write("
      data{
      int<lower=1> N;
      vector[N] x;
      vector[N] y;
      }
      
parameters{
    real alpha;
    real beta;
    real<lower=0> sigma;
}

model{
  alpha ~ normal(0,10);
  beta ~ normal(2,5);
  y ~ normal(alpha+beta*x, sigma);
}
generated quantities{
  real y_rep[N];
  for(n in 1:N) {
  y_rep[n] = normal_rng(x[n]*beta + alpha, sigma);
  }

}
   ",
"stan_model1.stan")

# compile the stan model:
stanc("stan_model1.stan")
stan_model1 <- "stan_model1.stan"

# do the regression with the stanmodel
fit <- stan(file = stan_model1, data = stan_data, warmup = 500, iter = 1000, chains = 4)
fit
posterior <- extract(fit)
plot(furan ~time, pch=20)
abline(mean(posterior$alpha),mean(posterior$beta), col = 6, lw = 2)

# traceplots:
plot(posterior$alpha, type = "l")
plot(posterior$beta, type="l")
plot(posterior$sigma, type="l")

# another possibility:
traceplot(fit)

# densities
par(mfrow=c(1,3))
plot(density(posterior$alpha), main = "alpha")
plot(density(posterior$beta), main = "beta")
plot(density(posterior$sigma), main = "sigma")

stan_dens(fit)
stan_hist(fit)

y_rep<- as.matrix(fit, pars = "y_rep")
dim(y_rep)
str(y_rep)
library(bayesplot)
ppc_dens_overlay(furan, y_rep[1:200,])
ppc_stat(y=furan, yrep=y_rep, stat="mean")
ppc_scatter_avg(y=furan, yrep=y_rep)
