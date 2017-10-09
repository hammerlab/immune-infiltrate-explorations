// modified to feed in priors
// attempting to put in vector rather than int
// update: now trying to find the available argument signatures for matrix

data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  matrix[3,3] alpha_mu;
  matrix[3,3] alpha_std;
}
parameters {
  matrix[3,3] alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(alpha_mu,alpha_std);    
  beta ~ normal(0,10);
  sigma ~ cauchy(0,5);
  for (n in 1:N)
    y[n] ~ normal(alpha[1,1] + beta * x[n], sigma);
}