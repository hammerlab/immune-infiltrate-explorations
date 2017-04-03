# use poisson distribution
# re-write to use matrix notation
data {
    // dimensions
    int<lower=1> N;  // N obs
    int<lower=1> G;  // N genes
    int<lower=1> S;  // N samples
    int<lower=0> C;  // N classes (e.g. B-cell, T-cell, B_Naive, CD5, CD45RO, etc)
                     // note: classes should be mutually exclusive. Each row here should sum to 1
   
    // data
    int<lower=1, upper=G> gene[N];    // gene id for each obs
    int<lower=1, upper=S> sample[N];  // sample id for each obs
    vector<lower=0, upper=1>[C] x[N]; // map each obs to each class (0:'- or ?', 1:'+')
    int<lower=0> y[N];                // count/tpm for each obs
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
    matrix<lower=0>[G, C] theta;   // loading factors for each gene, for each cell type
    vector[G] log_gene_base;       // constant intercept expression level for each gene, irrespective of cell type
}
model {
    for (i in 1:G)
        theta[i] ~ normal(0, 1);
    log_gene_base ~ normal(0, 1);
    for (s in 1:S) {
        vector[G] log_expected_rate;
        log_expected_rate = log_gene_base + log(theta*sample_x[s]);
        sample_y[s] ~ poisson_log(log_expected_rate);
    }
}
generated quantities {
    int y_rep[N];
    real log_lik[N];
    
    for (n in 1:N) {
        real log_expected_rate;
        log_expected_rate = log_gene_base[gene[n]] + log(theta[gene[n], ]*x[n]);
        y_rep[n] = poisson_log_rng(log_expected_rate);
        log_lik[n] = poisson_log_lpmf(y[n] | log_expected_rate);
    }
}