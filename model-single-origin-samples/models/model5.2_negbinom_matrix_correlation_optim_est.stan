## neg binom parameterization
## estimate correlation matrix among cell types
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
    
    // out of sample estimates, with unknown comp
    int<lower=1> N2;           // number of records in out of sample (UNK)
    int<lower=1> S2;           // number of samples in UNK set
    int<lower=1, upper=G> gene2[N2];    // gene id for UNK data (corresponding to IDs above)
    int<lower=1, upper=S2> sample2[N2]; // sample id for each UNK sample (separate from above)
    int<lower=0> y2[N2];       // data for UNK set 
}
transformed data {
    int sample_y[S, G];    // array (size SxG) of ints
    vector[C] sample_x[S]; // array (size S) of vectors[C]
    int sample2_y[S2, G];
    vector[C] theta_mu;
    
    for (n in 1:N) {
        sample_y[sample[n], gene[n]] = y[n];
        sample_x[sample[n]] = x[n];
    }
    for (n in 1:N2) {
        sample2_y[sample2[n], gene2[n]] = y2[n];
    }

    for (c in 1:C)
        theta_mu[c] = 0;
}
parameters {
    corr_matrix[C] Omega;        // degree of correlation among loading factors for each cell type
    vector<lower=0>[C] tau;      // scale for each cell type - multiplied (on diagonal) with Omega
    matrix<lower=0>[G, C] theta; // loading factors for each gene, for each cell type
    //vector[C] theta_mu;
    
    vector[G] log_gene_base;     // constant intercept expression level for each gene, irrespective of cell type
    vector<lower=0>[G] gene_phi; // overdispersion parameter per transcript (for now)
    
    simplex[C] sample2_x[S2];     // inferred sample2 compositions (simplex type enforces sum-to-one)
}
model {
    tau ~ cauchy(0, 2.5);
    Omega ~ lkj_corr(2);
    //theta_mu ~ normal(0, 1);
    for (g in 1:G) {
        theta[g] ~ multi_normal(theta_mu, quad_form_diag(Omega, tau));
    }
    log_gene_base ~ normal(0, 1);
    gene_phi ~ normal(0, 1);
    for (s in 1:S) {
        vector[G] log_expected_rate;
        log_expected_rate = log_gene_base + log(theta*sample_x[s]);
        sample_y[s] ~ neg_binomial_2_log(log_expected_rate, gene_phi);
    }
    for (s in 1:S2) {
        vector[G] log_expected_rate;
        log_expected_rate = log_gene_base + log(theta*sample2_x[s]);
        sample2_y[s] ~ neg_binomial_2_log(log_expected_rate, gene_phi);
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
