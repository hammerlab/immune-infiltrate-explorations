// neg binom parameterization
// estimate correlation matrix among cell types

// other bucket
// test first with other proportion of 0. i.e. with synthetic mixtures of known cell types
// then simulate other data:
// either pick other genes that we call tumor-related and give some values for those and spike in
// or spike in cell line or RCC sample (though that may have immune content)

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
    int<lower=1> nu;
    for (n in 1:N) {
        sample_y[sample[n], gene[n]] = y[n];
        sample_x[sample[n]] = x[n,];
    }
    for (n in 1:N2) {
        sample2_y[sample2[n], gene2[n]] = y2[n];
    }
    nu = 1;
}
parameters {
    cholesky_factor_corr[C] Omega_L;
    vector<lower=0>[C] Omega_sigma;
    //corr_matrix[C] Omega;        // degree of correlation among loading factors for each cell type
    //vector<lower=0>[C] tau;      // scale for each cell type - multiplied (on diagonal) with Omega
    
    matrix<lower=0>[G, C] theta; // loading factors for each gene, for each cell type
    vector[C] theta_mu;          // mean expression level for each cell type
    vector[M] theta_coefs_raw;
    vector[M] theta_coefs_per_gene[G];
    
    vector[G] log_gene_base;     // constant intercept expression level for each gene, irrespective of cell type
    // vector<lower=0>[G] gene_phi; // overdispersion parameter per transcript (for now)
    simplex[C] sample2_x[S2];     // inferred sample2 compositions (simplex type enforces sum-to-one)

    vector<lower=0, upper=1>[S2] unknown_prop; // proportion of each test sample that is of unknown cell type
    vector[G] other_log_contribution_per_gene[S2]; // for each test sample, per-transcript contribution of unknown cell type

    // hierarchical dispersion
    real log_global_phi_scale;      // single global-scale value that doesn't vary by sample, gene, or cell type.
    simplex[G] log_gene_phi;        // constant intercept overdispersion parameter per transcript. simplex solves identifiability problem
    vector[S] log_sample_scale;     // sample-level variance
    vector[S2] log_sample2_scale;   // sample2-level variance
    matrix[G, C] celltype_scale;    // cell type-level variance (varies by transcript)

}
transformed parameters {
    vector[M] theta_coefs[G];
    vector[G] corrected_phis[S];
    vector[G] corrected_phis2[S2];
    vector[G] ones;
    
    for (g in 1:G) 
        theta_coefs[g] = theta_coefs_raw + theta_coefs_per_gene[g];

    for (s in 1:S)
        corrected_phis[s] = log_global_phi_scale + log_gene_phi + log_sample_scale[s] + log(celltype_scale * sample_x[s]);
    corrected_phis = exp(corrected_phis);

    for (s in 1:S2)
        corrected_phis2[s] = log_global_phi_scale + log_gene_phi + log_sample2_scale[s] + log(celltype_scale * sample_x[s]);
    corrected_phis2 = exp(corrected_phis2);

    for(g in 1:G)
        ones[g] = 1; // initialize a vector of all ones (for dirichlet prior)
}
model {
    // estimate theta - gene-level expression per cell type, as a function of cell-surface expression proteins
    theta_mu ~ normal(0, 1);
    theta_coefs_raw ~ normal(0, 1);
    Omega_sigma ~ gamma(0.1, 0.1);
    Omega_L ~ lkj_corr_cholesky(nu);
    for (g in 1:G) {
        vector[C] theta_tmp; // temporary predictor for cell-gene-specific expression level
        theta_coefs_per_gene[g] ~ normal(0, 1);
        theta_tmp = theta_mu + cell_features*theta_coefs[g];
        theta[g] ~ multi_normal_cholesky(theta_tmp, diag_pre_multiply(Omega_sigma, Omega_L));
    }

    // estimate sample_y: observed expression for a sample (possibly a mixture)
    log_gene_base ~ normal(0, 1);
    log_global_phi_scale ~ cauchy(0,1); // fatter tails so that this can move more.
    log_gene_phi ~ dirichlet(ones);
    log_sample_scale ~ normal(0, 1);
    log_sample2_scale ~ normal(0, 1);
    for(g in 1:G)
        celltype_scale[g] ~ normal(0, 1);

    for (s in 1:S) {
        vector[G] log_expected_rate;
        log_expected_rate = log_gene_base + log(theta*sample_x[s]);
        sample_y[s] ~ neg_binomial_2_log(log_expected_rate, corrected_phis[s]);
    }
    
    // estimate sample2_y: observed expression for a sample of unknown composition
    unknown_prop ~ beta(5, 5); // not sure about this, maybe Beta(1,1) uniform?
    

    for (s in 1:S2) {
        vector[G] log_expected_rate;
        other_log_contribution_per_gene[s] ~ normal(0, 1);
        // TODO: use logmix() instead after `log_gene_base + ` in next line.
        log_expected_rate = log_gene_base + log(theta*sample2_x[s]) * (1 - unknown_prop[s]) + other_log_contribution_per_gene[s] * unknown_prop[s];
        sample2_y[s] ~ neg_binomial_2_log(log_expected_rate, corrected_phis2[s]);
    }
}
generated quantities {
    int y_rep[N];
    real log_lik[N];
    matrix[C,C] Omega;
    matrix[C,C] tau;
    Omega = multiply_lower_tri_self_transpose(Omega_L);
    tau = quad_form_diag(Omega_L, Omega_sigma);
    
    for (n in 1:N) {
        real log_expected_rate;
        log_expected_rate = log_gene_base[gene[n]] + log(theta[gene[n], ]*x[n]);
        y_rep[n] = neg_binomial_2_log_rng(log_expected_rate, corrected_phis[sample[n]][gene[n]]);
        log_lik[n] = neg_binomial_2_log_lpmf(y[n] | log_expected_rate, corrected_phis[sample[n]][gene[n]]);
    }
}