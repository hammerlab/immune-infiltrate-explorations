
# coding: utf-8

# In[1]:

import data
import models
import cache
import seaborn as sns
import numpy as np
import pandas as pd
import stanity
from stancache import stancache
import plotly.plotly as py
import plotly.graph_objs as go


# In[2]:

sns.set(context='talk')


# In[12]:

model_name = 'model5.3'
by = 'SubSet'
sample_n = 100
cache_only = False


# ## get data, as we did in earlier examples

# This will help in case we want to compare estimates for particular genes or samples

# In[4]:

sample_df = cache.cached(models.prep_sample_df, sample_n=sample_n)


# In[5]:

stan_data1 = models.prep_stan_data(sample_df, by=by, nu=1)
stan_data2 = models.prep_stan_data(sample_df, by=by, nu=2)
stan_data3 = models.prep_stan_data(sample_df, by=by, nu=3)
stan_data4 = models.prep_stan_data(sample_df, by=by, nu=4)
stan_data5 = models.prep_stan_data(sample_df, by=by, nu=5)
stan_data6 = models.prep_stan_data(sample_df, by=by, nu=6)


# In[6]:

model = models.get_model_file(model_name=model_name)


# ## get models from cache

# In[7]:

fit_nu1 = models.cached_stan_fit(file=model, data=stan_data1,
                              model_name=model_name, cache_only=cache_only)
fit_nu1_prefix = stancache.cached_stan_file(file=model, data=stan_data1, model_name=model_name, prefix_only=True)
fit_nu1_prefix


# In[8]:

fit_nu2 = models.cached_stan_fit(file=model, data=stan_data2, model_name=model_name, cache_only=cache_only)
fit_nu2_prefix = stancache.cached_stan_file(file=model, data=stan_data2, model_name=model_name, prefix_only=True)
fit_nu2_prefix


# In[9]:

fit_nu3 = models.cached_stan_fit(file=model, data=stan_data3, model_name=model_name, cache_only=cache_only)
fit_nu3_prefix = stancache.cached_stan_file(file=model, data=stan_data3, model_name=model_name, prefix_only=True)
fit_nu3_prefix


# In[10]:

fit_nu4 = models.cached_stan_fit(file=model, data=stan_data4, model_name=model_name, cache_only=cache_only)
fit_nu4_prefix = stancache.cached_stan_file(file=model, data=stan_data4, model_name=model_name, prefix_only=True)
fit_nu4_prefix


# In[11]:

fit_nu5 = models.cached_stan_fit(file=model, data=stan_data5, model_name=model_name, cache_only=cache_only)
fit_nu5_prefix = stancache.cached_stan_file(file=model, data=stan_data5, model_name=model_name, prefix_only=True)
fit_nu5_prefix


# In[13]:

fit_nu6 = models.cached_stan_fit(file=model, data=stan_data6, model_name=model_name, cache_only=cache_only)
fit_nu6_prefix = stancache.cached_stan_file(file=model, data=stan_data6, model_name=model_name, prefix_only=True)
fit_nu6_prefix


# ## compute loo-psis for each model

# In[14]:

loo_nu1 = cache.cached(stanity.psisloo,
                        log_likelihood=fit_nu1.extract('log_lik')['log_lik'],
                    cache_filename='{}.loo.pkl'.format(fit_nu1_prefix))
loo_nu1.print_summary()


# In[15]:

loo_nu2 = cache.cached(stanity.psisloo,
                      log_likelihood=fit_nu2.extract('log_lik')['log_lik'],
                      cache_filename='{}.loo.pkl'.format(fit_nu2_prefix))
loo_nu2.print_summary()


# In[16]:

loo_nu3 = cache.cached(stanity.psisloo,
                      log_likelihood=fit_nu3.extract('log_lik')['log_lik'],
                      cache_filename='{}.loo.pkl'.format(fit_nu3_prefix))
loo_nu3.print_summary()


# In[17]:

loo_nu4 = cache.cached(stanity.psisloo,
                      log_likelihood=fit_nu4.extract('log_lik')['log_lik'],
                      cache_filename='{}.loo.pkl'.format(fit_nu4_prefix))
loo_nu4.print_summary()


# In[18]:

loo_nu5 = cache.cached(stanity.psisloo,
                      log_likelihood=fit_nu5.extract('log_lik')['log_lik'],
                      cache_filename='{}.loo.pkl'.format(fit_nu5_prefix))
loo_nu5.print_summary()


# In[19]:

loo_nu6 = cache.cached(stanity.psisloo,
                      log_likelihood=fit_nu6.extract('log_lik')['log_lik'],
                      cache_filename='{}.loo.pkl'.format(fit_nu6_prefix))
loo_nu6.print_summary()


# ## compare psis-loo for nu==1 vs nu==3

# In[20]:

stanity.loo_compare(loo_nu1, loo_nu3)


# It seems that this model works better with higher values of `nu`.

# ## compare psis-loo for nu==2 vs nu==3

# In[21]:

stanity.loo_compare(loo_nu2, loo_nu3)


# ## compare psis-loo for nu==3 vs nu==4

# In[22]:

stanity.loo_compare(loo_nu3, loo_nu4)


# ## compare psis-loo for nu==4 vs nu==5

# In[23]:

stanity.loo_compare(loo_nu4, loo_nu5)


# ## compare psis-loo for nu==5 vs nu==6

# In[24]:

stanity.loo_compare(loo_nu5, loo_nu6)


# ## Summarize fit across different values of `nu`

# In[25]:

loo_results = {1: loo_nu1, 2: loo_nu2, 3: loo_nu3, 4: loo_nu4, 5: loo_nu5, 6: loo_nu6}


# In[26]:

loodf = list()
for (nu, loores) in loo_results.items():
    thisdf = loores.pointwise.reset_index()
    thisdf['nu'] = nu
    loodf.append(thisdf)
loodf = pd.concat(loodf)
loodf.head()


# In[27]:

loodf.sort_values(['index','nu'], inplace=True)
loodf['first_diff'] = loodf.groupby('index')['pointwise_elpd'].diff()
loodf.fillna(0, inplace=True)
loodf['cum_diff'] = loodf.groupby('index')['first_diff'].cumsum()
loodf.head()


# In[28]:

loodf2 = pd.merge(sample_df, loodf, on='index')


# ## plot ELPD by `nu`

# In[29]:

sns.regplot(data=loodf2, x='nu', y='first_diff', x_estimator=np.sum)


# In[30]:

sns.regplot(data=loodf2, x='nu', y='cum_diff', x_estimator=np.sum)


# In[31]:

by_nu = loodf2.groupby(['nu'])['cum_diff'].sum().reset_index()
trace = go.Scatter(
    x=by_nu['nu'],
    y=by_nu['cum_diff'],
    mode='markers+lines',
    line=dict(
        shape='spline'
    ),
)
py.iplot([trace])


# ## plot ELPD by observation over different values of `nu`

# In[32]:

from matplotlib import pyplot as plt
for i in np.arange(max(loodf.index)):
    plt.plot(loodf.loc[loodf['index'] == i, 'nu'],
             loodf.loc[loodf['index'] == i, 'first_diff'],
             alpha=0.2)


# ## plot `elpd` by cell_type

# In[33]:

by_subset = loodf2.groupby(['cell_type','nu']).agg({'first_diff': np.sum, 'cum_diff': np.sum}).reset_index()
plot_data = list()
for (name, sample_df) in by_subset.groupby('cell_type'):
    trace = go.Scatter(
        x=sample_df['nu'],
        y=sample_df['cum_diff'],
        mode='markers+lines',
        name=name,
        line=dict(
            shape='spline'
        ),
    )
    plot_data.append(trace)
py.iplot(plot_data, filename='immune-inf/elpd-by-cell_type')


# ## plot elpd by `SubSet`
# 
# Is there any pattern in how errors are distributed by SubSet?

# In[34]:

by_subset = loodf2.groupby(['SubSet','nu']).agg({'first_diff': np.sum, 'cum_diff': np.sum}).reset_index()
plot_data = list()
for (name, sample_df) in by_subset.groupby('SubSet'):
    trace = go.Scatter(
        x=sample_df['nu'],
        y=sample_df['cum_diff'],
        mode='markers+lines',
        name=name,
        line=dict(
            shape='spline'
        ),
    )
    plot_data.append(trace)
py.iplot(plot_data, filename='immune-inf/elpd-by-subset')


# ## look at elpd over `nu`, by sample_id

# In[35]:

by_sample = loodf2.groupby(['sample_id', 'SubSet', 'nu'])['first_diff'].sum().reset_index()
plot_data = list()
for (sample_id, sample_df) in by_sample.groupby('sample_id'):
    trace = go.Scatter(
        x=sample_df['nu'],
        y=sample_df['first_diff'],
        mode='lines',
        name=sample_df['SubSet'].values[0],
        text=sample_id,
        hoverinfo='text+name',
        line=dict(
            shape='linear'
        ),
        legendgroup=sample_df['SubSet'].values[0],
        opacity=0.2,
    )
    plot_data.append(trace)


# In[36]:

py.iplot(plot_data)


# In[ ]:



