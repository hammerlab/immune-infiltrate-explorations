## neg binom parameterization
## estimate correlation matrix among cell types
data {
    // dimensions
    int<lower=1> N;  // N obs
    int<lower=1> G;  // N genes
    int<lower=1> S;  // N samples
    int<lower=0> C;  // N classes (e.g. B-cell, T-cell, B_Naive, CD5, CD45RO, etc)
                     //     note: classes should be mutually exclusive. Each row here should sum to 1
    int<lower=0> M; // number of cell-level predictors 
   
    // data for each gene*sample
    int<lower=1, upper=G> gene[N];    // gene id for each obs
    int<lower=1, upper=S> sample[N];  // sample id for each obs
    vector<lower=0, upper=1>[C] x[N]; // map each obs to each class (0:'- or ?', 1:'+')
    int<lower=0> y[N];                // count/tpm for each obs
    
    // group-level predictors for each class C
    matrix[C, M] cell_features; 
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
    //corr_matrix[C] Omega;        // degree of correlation among loading factors for each cell type
    vector<lower=0>[C] tau;        // scale for each cell type - multiplied (on diagonal) with Omega
    matrix[M, G] theta_coefs_per_gene;
    
    vector[G] log_gene_base;     // constant intercept expression level for each gene, irrespective of cell type
    vector<lower=0>[G] gene_phi; // overdispersion parameter per transcript (for now)
}
transformed parameters {
    matrix<lower=0>[G, C] theta;
    { 
        matrix[C, G] tmp_theta;
        for (c in 1:C)
            tmp_theta[c] = exp(cell_features[c]*theta_coefs_per_gene); // better to do on log scale?
        theta = tmp_theta' + (diag_pre_multiply(tau, L_Omega) * z)';
    }
}
model {
    // priors on components of theta: relative expression per cell type per gene transcript
    tau ~ cauchy(0, 2.5);
    to_vector(z) ~ normal(0, 1);
    L_Omega ~ lkj_corr_cholesky(2);
    to_vector(theta_coefs_per_gene) ~ normal(0, 1);
    
    // estimate sample_y: obseved expression for a sample (possibly a mixture)
    log_gene_base ~ normal(0, 1);
    gene_phi ~ normal(0, 1);
    for (s in 1:S) {
        vector[G] log_expected_rate;
        log_expected_rate = log_gene_base + log(theta*sample_x[s]);
        sample_y[s] ~ neg_binomial_2_log(log_expected_rate, gene_phi);
    }
}
generated quantities {
    int y_rep[N];
    real log_lik[N];
    
    for (n in 1:N) {
        real log_expected_rate;
        log_expected_rate = log_gene_base[gene[n]] + log(theta[gene[n], ]*x[n]);
        y_rep[n] = neg_binomial_2_log_rng(log_expected_rate, gene_phi[gene[n]]);
        log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_expected_rate, gene_phi[gene[n]]);
    }
}