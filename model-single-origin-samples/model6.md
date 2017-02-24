## Model 6

#### Likelihood

The input data is a matrix of gene counts $\boldsymbol Y$, a $(S \times G)$ matrix where entry $Y_{s,g}$ corresponds to the mean number (tpm) of transcripts of $g$ in sample $s$. 

`model6` is a generalized linear model. (N.B. these immune infiltrate explorations are not restricted to linear models and we explore nonlinear ones in the future.) Per the GLM requirement that the expectation of our outcome variable $y$ is a nonlinear function of the linear predictor,  

$$ E(y | X, \beta ) = g^{-1}(X\beta) $$

We choose the log of a negative binomial as $g$: 

$$ \vec{y_s} = NegBinLog( \vec{\mu} , \vec{\phi} )  $$

where $\mu$ is the log of the mean expression and $\phi$ the dispersion parameter 

$ \vec{\mu_{G \times 1}} = \vec{\tilde{\mu}_{G \times 1}} + log(\boldsymbol{\beta}_{G \times C} * \vec{x_{s_{C \times 1}}}) $ 


#### Prior on coefficient, $\beta$

$ \vec{\beta_g} = MultiNormal(\vec{u_g}, \boldsymbol \Sigma)$ 

where 

$ \vec{u_{g_{C \times 1}}} = \vec{p_{C \times 1}} + \boldsymbol{F}_{C \times M} * 
( \vec {b_{M \times 1}} + \vec{\kappa_{g_{M \times 1}}} ) $

and the covariance matrix $\Omega$ is decomposed into the diagonal matrix ("scaling factor") $\tau$ and correlation matrix $\Omega$

$$ \Sigma = \tau \space \Omega \space \tau $$




#### Hyper priors 

$ \tau \sim Cauchy(0, 2.5) $ Weakly informative, with small scale, and positive

$ \Omega \sim LKJcorr(\nu) $, here we chose $\nu = 2$

$ p \sim Normal(0,1) $

$ b \sim Normal(0, 1) $

$  \kappa_g \sim Normal(0,1) $

$ \tilde{\alpha} \sim Normal(0,1) $

$ \phi_{(G \times 1)} \sim Normal(0,1) $


#### Translation	

| J code  | E notation  | Description  |
|---|---|---|
| `theta`  |   $\beta$ |  loading factors for each gene x cell type |
| `theta_mu`   | $p$  |  mean expr. level for each cell type |
| `theta_coefs`  |  $b$ | coefficient per feature |
| `theta_coefs_per_gene`  |  $\kappa$ | Feature's influence on gene  |
| `theta_tmp` |  $u$ |   |
|  `log_expected_rate` | $\alpha$  |   |
|  `log_gene_base` | $\tilde{\alpha}$  |   |
| `cell_features` | $F$ | Cell-type features (for now just cell surface markers) |
| `gene_phi ` | $\phi$ | |

------

## Proposed extensions for `upenn-tex`

### 1. Incorporate a "gene signature" into model

E.g. in `upenn-tex`, we have up- and down-regulated gene sets that characterize the transcriptional state of exhausted T cells. We may easily incorporate these gene sets into our model by extending the matrix $\boldsymbol{F}$, which currently only holds cell surface marker features. 


#### 2. Incorporate TF information into model

We expect that transcripts regulated by the same transcription factors will be tightly co-expressed. We have two ways of incorporating this prior information into our model. 

1) One method of using $Inferelator$ is to first transform gene count data into *gene cluster* counts, where the level of each gene cluster is the arithmetic mean of the genes associated with it. This is helpful to them for reducing matrix sparsity. If we have information on a given cell type's gene regulatory network, a first experiment to assess the 

2) In `model6` we estimate the covariance matrix among expression profiles of different cell types. We may similarly estimate the covariance matrix among genes *within* each cell type, using the GRNs as prior knowledge.

### Directions 

1. How well can we do on inferring the presence of one well-characterized immune cell type in a known mixture?
	* In an unknown mixture? 
2. How well can we do on a less well-characterized immune cell type in a {known, unknown} mixture? 
	* I.e. obscure or remove some of the features from (1) and compare performance



