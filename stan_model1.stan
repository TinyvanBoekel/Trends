
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
  beta ~ normal(0,5);
  y ~ normal(alpha+beta*x, sigma);
}
generated quantities{
  real y_rep[N];
  for(n in 1:N) {
  y_rep[n] = normal_rng(x[n]*beta + alpha, sigma);
  }

}
   
