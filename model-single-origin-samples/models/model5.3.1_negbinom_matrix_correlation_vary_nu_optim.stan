## neg binom parameterization
## estimate correlation matrix among cell types
functions {
  /* 
    ## per : https://groups.google.com/forum/?fromgroups#!searchin/stan-users/neg_binomial_2_log_rng/stan-users/rgpHgBtFS3k/A7GqglKOAQAJ
    ## neg_binomial_2_rng(mu,psi) as suggested by Ben: 
    real gamma_temp; 
    if (is_inf(psi)) 
      gamma_temp <- mu; 
    else 
      gamma_temp <- gamma_rng(psi, psi / mu); 
    if (gamma_temp < pow(2.0, 30.0)) 
      result <- poisson_rng(gamma_temp); 
    else 
      result <- normal_rng(gamma_temp, sqrt(gamma_temp)); 
  */
  real alt_neg_binom_2_rng(real mu, real psi) {
    real gamma_temp; 
    real result;
    
    if (is_nan(psi)) 
      gamma_temp = mu;
    else if (is_inf(psi)) 
      gamma_temp = mu; 
    else if (is_nan(psi / mu))
      gamma_temp = mu;
    else if (is_inf(psi / mu))
      gamma_temp = mu;
    else if (psi / mu < 0)
      gamma_temp = mu;
    else 
      gamma_temp = gamma_rng(psi, psi / mu);
      
    if (is_nan(gamma_temp) || is_inf(gamma_temp) || gamma_temp < 0)
      result = -1;
    else if (gamma_temp < pow(2.0, 30.0)) 
      result = poisson_rng(gamma_temp);
    else 
      result = normal_rng(gamma_temp, sqrt(gamma_temp)); 
    
    if (result > 2^31) {
      result = 2^31;
    }
    
    return floor(result);
  }
}
data {
    // dimensions
    int<lower=1> N;  // N obs
    int<lower=1> G;  // N genes
    int<lower=1> S;  // N samples
    int<lower=0> C;  // N classes (e.g. B-cell, T-cell, B_Naive, CD5, CD45RO, etc)
                     //     note: classes should be mutually exclusive. Each row here should sum to 1
    // int<lower=0> M; // number of cell-level predictors 
   
    // data for each gene*sample
    int<lower=1, upper=G> gene[N];    // gene id for each obs
    int<lower=1, upper=S> sample[N];  // sample id for each obs
    vector<lower=0, upper=1>[C] x[N]; // map each obs to each class (0:'- or ?', 1:'+')
    int<lower=0> y[N];                // count/tpm for each obs
    
    int<lower=1> nu;  // hyper-parameter for lkj_corr prior in Omega.
    
}
transformed data {
    int sample_y[S, G];    // array (size SxG) of ints
    vector[C] sample_x[S]; // array (size S) of vectors[C]
    for (n in 1:N) {
        sample_y[sample[n], gene[n]] = y[n];
        sample_x[sample[n]] = x[n,];
    }
}
parameters {
    cholesky_factor_corr[C] L_Omega;
    matrix[C, G] z;
    vector<lower=0>[C] tau;        // scale for each cell type - multiplied (on diagonal) with Omega
    
    vector[C] theta_mu;
    vector[G] log_gene_base;       // constant intercept expression level for each gene
                                   // irrespective of cell type
    vector<lower=0>[G] gene_phi;   // overdispersion parameter per transcript (for now)
}
transformed parameters {
    matrix<lower=0>[G, C] theta;
    {
        matrix[G, C] log_theta;
        log_theta = rep_matrix(theta_mu', G) + (diag_pre_multiply(tau, L_Omega) * z)';
        theta = exp(log_theta);
    }
}
model {
    tau ~ cauchy(0, 2.5);
    to_vector(z) ~ normal(0, 1);
    L_Omega ~ lkj_corr_cholesky(nu);
    theta_mu ~ normal(0, 1);
    log_gene_base ~ normal(0, 1);
    gene_phi ~ normal(0, 1);
    for (s in 1:S) {
        vector[G] log_expected_rate;
        log_expected_rate = log_gene_base + log(theta*sample_x[s]);
        sample_y[s] ~ neg_binomial_2_log(log_expected_rate, gene_phi);
    }
}
generated quantities {
    real y_rep[N];
    real log_lik[N];
    
    for (n in 1:N) {
        real log_expected_rate;
        log_expected_rate = log_gene_base[gene[n]] + log(theta[gene[n], ]*x[n]);
        y_rep[n] = alt_neg_binom_2_rng(exp(log_expected_rate), gene_phi[gene[n]]);
        log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_expected_rate, gene_phi[gene[n]]);
    }
}
