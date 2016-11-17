
# coding: utf-8

# In[1]:

import data
import models
import cache
import seaborn as sns
import numpy as np
import pandas as pd
import patsy
from matplotlib import pyplot as plt


# In[2]:

sns.set(context='talk')


# In[3]:

model_name = 'model5'
by = 'SubSet'
sample_n = 1000


# ## load files for all cell types

# In[4]:

df = cache.cached(data.prep_annotated_data)


# In[5]:

assert all(pd.notnull(df['log1p_tpm_rescaled']))


# In[6]:

print(df.columns)

#apply(lambda x: x.startswith('C'))


# ## sample genes for analysis

# In[7]:

sample_df = cache.cached(models.prep_sample_df, sample_n=sample_n)


# ## fit model

# In[8]:

stan_data = models.prep_stan_data(sample_df, by=by)


# In[9]:

model_file = models.get_model_file(model_name=model_name)
print(cache._read_file(model_file))


# In[10]:

model_fit = models.cached_stan_fit(file=model_file, data=stan_data, model_name=model_name)


# ## check convergence (superficially)

# In[11]:

models.print_stan_summary(model_fit, pars='lp__')


# In[12]:

models.plot_stan_summary(model_fit, pars='theta', metric='Rhat')


# ## summarize posterior draws of theta by gene

# In[13]:

# meta-data used for plotting functions below
# so that the following code is invariant to the model run
colnames = list(stan_data['x'].columns)
sort_by = colnames[0]
print(sort_by)


# In[14]:

theta_ldf = models.prep_theta_summary(model_fit, sample_df=sample_df, colnames=colnames, expose_group=sort_by)


# In[15]:

## show theta estimates for first 50 genes, by `sort-by`
g = sns.boxplot(data=theta_ldf.loc[theta_ldf['mean_value_rank_{}'.format(sort_by)] <= 50,:]                 .sort_values('mean_value_rank_{}'.format(sort_by)),
            y='new_gene_cat',
            x='value',
            hue='variable', 
            fliersize=0, width=2, linewidth=0.2)


# In[16]:

## zoom in on the highest-ranked genes by `sort-by` difference from average 
## across all cell types
g = sns.boxplot(data=theta_ldf.loc[theta_ldf['mean_abs_diff_rank_{}'.format(sort_by)] <= 10,:]                 .sort_values('mean_diff_rank_{}'.format(sort_by)),
            y='new_gene_cat',
            x='value',
            hue='variable', 
            fliersize=0, linewidth=0.2)


# ## posterior-predictive checking for selected genes

# In[17]:

# identify top_genes by name
top_genes = theta_ldf.loc[theta_ldf['mean_abs_diff_rank_{}'.format(sort_by)] <= 10,:]                 .drop_duplicates(subset='new_gene_cat')['new_gene_cat'].values
print(top_genes)


# In[18]:

# get yrep draws
yrep_df = models.prep_yrep_summary(model_fit, sample_df=sample_df, filter_genes=top_genes)


# In[19]:

models.plot_posterior_predictive_checks(model_fit=model_fit, plot_genes=top_genes, sample_df=sample_df,
                                        yrep_df=yrep_df, n_genes=1)


# In[20]:

models.plot_posterior_predictive_checks(model_fit=model_fit, plot_genes=top_genes[1:], sample_df=sample_df,
                                        yrep_df=yrep_df, n_genes=1)


# In[21]:

models.plot_posterior_predictive_checks(model_fit=model_fit, plot_genes=top_genes[2:], sample_df=sample_df,
                                        yrep_df=yrep_df, n_genes=1)


# ## summarize posterior draws for `theta_mu`

# In[22]:

mu_ldf = models.prep_theta_mu_summary(stan_fit=model_fit, stan_data=stan_data, par='theta_mu')


# In[23]:

a = sns.boxplot(data=mu_ldf, x='variable', y='value')
plt.setp(a.get_xticklabels(), rotation='vertical')


# ## summarize posterior draws for `Omega`

# In[24]:

omega_summary = models.prep_omega_summary(stan_fit=model_fit, stan_data=stan_data, par='Omega', gene_id=by)


# In[25]:

with sns.plotting_context('paper'):
    sns.heatmap(omega_summary.loc[:, list(stan_data['x'].columns)])


# In[26]:



