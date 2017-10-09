// modified to feed in priors
// attempting to put in vector rather than int

data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  real alpha_mu[3];
  real alpha_std[3];
}
parameters {
  real alpha[3];
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(alpha_mu,alpha_std);    
  beta ~ normal(0,10);
  sigma ~ cauchy(0,5);
  for (n in 1:N)
    y[n] ~ normal(alpha[1] + beta * x[n], sigma);
}