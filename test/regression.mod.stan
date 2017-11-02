// modified to feed in priors

data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  int alpha_mu;
  int alpha_std;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(alpha_mu,alpha_std);    
  beta ~ normal(0,10);
  sigma ~ cauchy(0,5);
  for (n in 1:N)
    y[n] ~ normal(alpha + beta * x[n], sigma);
}