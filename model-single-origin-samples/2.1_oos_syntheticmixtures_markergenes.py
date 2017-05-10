
# coding: utf-8
# run this as: python -i scriptname.py | tee logfileout.txt

# In[1]:

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import data
import models
import cache
import pandas as pd
import matplotlib.pyplot as plt
#get_ipython().magic('matplotlib inline')
import seaborn as sns
import pystan
from time import time
from datetime import timedelta
import pickle
import dill

# In[2]:

bindea = pd.read_csv('../curated_data/pure_samples/bindea/Bindea Immunome - Bindea paper.csv')
bindea.head()


# In[3]:

marker_genes = bindea.Symbol.unique()
len(marker_genes)


# In[7]:

model_name = 'model6.2'
by = 'SubSet'
total_num_genes = 800 # number we want to have total


# to modify cache filename: if using `models.cached_stan_fit()`, change `model_name` parameter; if using `cache.cached()` directly, can change `file_prefix` parameter.

# In[16]:

sample_df = cache.cached(models.prep_sample_df,
                         sample_n=None)


# In[17]:

sample_df.SubSet.unique()


# In[18]:

# get all genes in sample_df
all_genes = sample_df.gene_name.unique()
len(all_genes)


# In[34]:

# how many marker genes are not in sample_df?
len(set(marker_genes) - set(all_genes))


# In[38]:

marker_genes_adjusted = list(set(marker_genes) - (set(marker_genes) - set(all_genes)))
len(marker_genes_adjusted)


# In[39]:

genes_to_sample_from = list(set(all_genes) - set(marker_genes_adjusted))
len(genes_to_sample_from)


# In[40]:

print('desired size:', total_num_genes-len(marker_genes_adjusted))
sampled_genes = np.random.choice(genes_to_sample_from, size=total_num_genes-len(marker_genes_adjusted))
len(sampled_genes)


# In[41]:

np.concatenate([marker_genes_adjusted, sampled_genes]).shape


# In[42]:

# following models.prep_sample_df() source:
filtered_sample_df = sample_df[sample_df['gene_name'].isin(np.concatenate([marker_genes_adjusted, sampled_genes]))].copy()
filtered_sample_df.sort_values(['gene_id','sample_id'], inplace=True)
filtered_sample_df['new_gene_cat'] = filtered_sample_df['gene_name'].astype('category')
filtered_sample_df['new_gene_id'] = filtered_sample_df['new_gene_cat'].cat.codes+1
filtered_sample_df['new_sample_cat'] = filtered_sample_df['sample_id'].astype('category')
filtered_sample_df['new_sample_id'] = filtered_sample_df['new_sample_cat'].cat.codes+1
filtered_sample_df.reset_index(inplace=True)
#assert len(filtered_sample_df.gene_name.unique()) == total_num_genes
# this fails but let's continue anyway

# In[43]:

sample_df.shape, filtered_sample_df.shape


# In[44]:

len(filtered_sample_df.gene_name.unique())
print('filtered sample df number of genes:', len(filtered_sample_df.gene_name.unique()))



# good enough!

# In[62]:

# make training_df = sample_df, test_df = some synthetic mixtures
# first let's try making these mixtures manually

"""
notes from before about how we want to do this:
take raw lines, add them straight up with weights
naive B cells, memory B cells
(CD8, CD4)
Tregs vs naive B cells

first level of noise: gaussian noise on every transcript after the sum
so take mix1 + np.random.normal(0,1,len(mix1))
doesn't feel like enough noise tbh

second level of noise: vary the mixture weights for each transcript individually
weights_noisy = np.reshape(weights*len(mix1), reference[simple_cols].shape)
weights_noisy += np.abs(np.random.normal(0,0.1, weights_noisy.shape))
weights_noisy = weights_noisy.clip(0,1)
mix3=(weights_noisy * reference[simple_cols]).sum(axis=1)
mix3.head()

"""

training_df = filtered_sample_df


# let's try running this thru stan now

# In[56]:

model_file = models.get_model_file(model_name=model_name)
#print(cache._read_file(model_file))


# In[63]:

stan_data = models.prep_stan_data(sample_df=training_df, test_df=None, by=by)
stan_data


# In[64]:

#relevant_sample_ids

def get_sample_ids_by_subset(training_df):
    return {subset: 
            sample_df[sample_df['SubSet'] == subset].new_sample_id.unique() for subset in sample_df.SubSet.unique()}
relevant_sample_ids = get_sample_ids_by_subset(training_df)
assert all([len(relevant_sample_ids[i]) > 0 for i in relevant_sample_ids])
relevant_sample_ids


# In[65]:

# wrap above manual work into a function

def mix_cell_lines(sample_df, xdata, subsets, weights, sample_ids=None, new_sample_id=10001):
    """
    e.g. xdata=stan_data['x'], subsets=['B_Naive', 'B_Memory'], weights=[.5, .5], sample_ids=None
    if sample_ids are None, the first sample of each subset is used
    """
    
    assert len(weights) == len(subsets)
    if not sample_ids:
        sample_ids = [relevant_sample_ids[subset][0] for subset in subsets]
    
    weights = np.array(weights)
    weights = weights / np.sum(weights) # normalize
    
    x2_data = pd.DataFrame(np.zeros((1, xdata.shape[1])), columns=xdata.columns)
    
    transformed_lines = []
    for subset, weight, sample_id in zip(subsets, weights, sample_ids):
        transformed = sample_df[sample_df['new_sample_id'] == sample_id].copy()
        transformed.loc[:,'est_counts'] *= weight
        transformed_lines.append(transformed)
        x2_data['SubSet[%s]' % subset] = weight
    
    mixed_sample = pd.concat(transformed_lines).groupby(         ['gene_name', 'new_gene_id'])['est_counts']        .sum().reset_index()
        
    mixed_sample['sample_id'] = new_sample_id
    
    return mixed_sample, x2_data
    
    
mix_cell_lines(sample_df=training_df, xdata=stan_data['x'], subsets=['B_Naive', 'B_Memory'], weights=[.5, .5], sample_ids=None)


# In[60]:

relevant_sample_ids


# In[66]:

# here are the mixtures we want
mix1, mix1_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'B_Memory'],
                              weights=[.5, .5],
                              sample_ids=[7, 4],
                              new_sample_id = 10001)

mix2, mix2_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'B_Memory'],
                              weights=[.5, .5],
                              sample_ids=[21, 29],
                              new_sample_id = 10002)

mix3, mix3_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'B_Memory'],
                              weights=[.25, .75],
                              sample_ids=[7, 4],
                              new_sample_id = 10003)

mix4, mix4_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'B_Memory'],
                              weights=[.25, .75],
                              sample_ids=[21, 29],
                              new_sample_id = 10004)

mix5, mix5_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'B_Memory'],
                              weights=[.75, .25],
                              sample_ids=[7, 4],
                              new_sample_id = 10005)

mix6, mix6_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'B_Memory'],
                              weights=[.75, .25],
                              sample_ids=[21, 29],
                              new_sample_id = 10006)

# tregs vs naive B cells

mix7, mix7_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'CD4_Treg'],
                              weights=[.5, .5],
                              sample_ids=[7, 18],
                              new_sample_id = 10007)

mix8, mix8_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'CD4_Treg'],
                              weights=[.25, .75],
                              sample_ids=[7, 18],
                              new_sample_id = 10008)

mix9, mix9_x = mix_cell_lines(sample_df=training_df,
                              xdata=stan_data['x'],
                              subsets=['B_Naive', 'CD4_Treg'],
                              weights=[.5, .5],
                              sample_ids=[21, 24],
                              new_sample_id = 10009)

mix10, mix10_x = mix_cell_lines(sample_df=training_df,
                                xdata=stan_data['x'],
                              subsets=['B_Naive', 'CD4_Treg'],
                              weights=[.25, .75],
                              sample_ids=[21, 24],
                              new_sample_id = 10010)


# In[67]:

# make a test_df and x2_data with all of them
test_df = pd.concat([mix1, mix2, mix3, mix4, mix5, mix6, mix7, mix8, mix9, mix10])
test_df['gene_id'] = test_df['new_gene_id']
x2_data = pd.concat([mix1_x,mix2_x,mix3_x,mix4_x,mix5_x,mix6_x,mix7_x,mix8_x,mix9_x,mix10_x])

#for dat in [small_training_df, test_df]:
for dat in [training_df, test_df]:
    dat.sort_values(['gene_id','sample_id'], inplace=True)
    dat['new_sample_cat'] = dat['sample_id'].astype('category')
    dat['new_sample_id'] = dat['new_sample_cat'].cat.codes+1


# In[68]:

test_df.head()


# In[69]:

x2_data.head()


# In[70]:

training_df.head()


# In[71]:

test_data = {'N2': len(test_df.index),
             'S2': len(test_df.new_sample_id.unique()),
             'gene2': test_df.new_gene_id.values,
             'sample2': test_df.new_sample_id.values,
             'y2': test_df.est_counts.astype(int).values,
             'x2': x2_data, ## for easy access later
             }
test_data


# In[72]:

stan_data.update(test_data)


# In[81]:

stan_data['G']


# In[74]:

model_name_for_cache = model_name + '.markergenes'
#model_fit_test = models.cached_stan_fit(file=model_file,
#                                   data=stan_data,
#                                   model_name=model_name_for_cache,
#                                   iter=5,
#                                   chains=1,
#                                   warmup=0,
#                                   refresh=25,
#                                   sample_file='logs/test_sampling_log.txt',
#                                   #diagnostic_file='logs/test_diagnostic_log.txt' # NotImplementedError
#                                  )


# In[79]:

#with open('2.1_extract_cpp_code.cpp', 'w') as w:
#    w.write(model_fit_test.stanmodel.get_cppcode())


# In[ ]:




# the sampling log was 6.7MB for 5 iterations!!! too much.
# 
# the real thing:

# In[ ]:

# model_fit = models.cached_stan_fit(file=model_file,
#                                    data=stan_data,
#                                    model_name=model_name_for_cache,
#                                    refresh=25,
#                                    sample_file='logs/markergenes_sampling_log.txt', # going to be huge
#                                    #iter=10000,
#                                    #chains=4,
#                                    #warmup=2500
#                                   )



# TODO: put in our manual pickling

print('Compile')
# stan_model = cached(func=pystan.StanModel,
#                      cache_filename=model_cachefile,
#                      model_code=model_code,
#                      cache_dir=cache_dir,
#                      model_name=model_name,
#                      cache_only=cache_only,
#                      force=force)

from stancache.stancache import _sanitize_model_name, _get_model_code
model_name_for_cache = _sanitize_model_name(model_name_for_cache)
model_code = _get_model_code(model_code=None, file=model_file)

start = time()
stan_model = pystan.StanModel(model_code=model_code,
                  model_name=model_name_for_cache
                )

end = time()
try:
  elapsed = str(timedelta(seconds=end-start))
  print('{}: Execution completed ({} elapsed)'.format('pystan.StanModel', elapsed))
  pickle.dump(stan_model, open('logs/stanmodel_markergenes.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)     
except:
  print('unexecpted error: {}'.format(sys.exc_info()[0]))

print('Get posterior draws from model')
start = time()
# fit = cached(func=stan_model.sampling,
#               cache_filename=fit_cachefile,
#               cache_dir=cache_dir,
#               force=force,
#               cache_only=cache_only,
#               **kwargs)

model_fit = stan_model.sampling(data=stan_data,
                    refresh=25,
                    sample_file='logs/markergenes_sampling_log.txt',
                    #iter=5,
                    #chains=1,
                    #warmup=0
)

end = time()
try:
  elapsed = str(timedelta(seconds=end-start))
  print('{}: Execution completed ({} elapsed)'.format('stan_model.sampling', elapsed))
  pickle.dump(model_fit, open('logs/modelfit_markergenes.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)     
except:
  print('unexecpted error: {}'.format(sys.exc_info()[0]))

print("Analysis...")

# ## inspect posterior predictions

# First, we note that the model didn't converge very well. We should probably diagnose fit of the model better.

# In[61]:

with open('logs/markergenes_lp.stansummary', 'w') as w:
  w.write(models.filter_stan_summary(model_fit, pars='lp__').to_string())


# However, the estimates for our sample2_x, which indicates mixing proportions of different types of cells, are somewhat better. Not great, but better.

# In[62]:

with open('logs/markergenes_sample2x.stansummary', 'w') as w:
  w.write(models.filter_stan_summary(model_fit, pars='sample2_x').to_string())


# These values `[0, x]` refer to the `[sample_id, subset_id]`. We only have one test sample, so this is fairly easy to interpret. 
# 
# Our subset_ids are indexed by the columns of the `stan_data['x']` object we passed in. 

# In[63]:

colnames = list(stan_data['x'].columns)
colnames


# We can re-use the code we wrote for `extract_theta_summary` to extract these values, since the data structure of `sample2_x` & of theta are similar. 

# In[64]:

inferred_type = models.extract_theta_summary(colnames=colnames, par='sample2_x', stan_fit=model_fit, gene_id='sample_id')


# In[65]:

inferred_type


# We then reshape & format our desired "subset" variable. 

# In[66]:

import re
df = pd.melt(inferred_type, id_vars=['iter','sample_id'], value_name='estimate', var_name='variable')
df['SubSet'] = df['variable'].apply(lambda x: re.sub(string=x, pattern='(.*)\[(.*)\]', repl='\\2'))


# In[67]:

df.head()


# In[68]:

df.shape


# In[69]:

df.groupby('sample_id').iter.count()


# In[70]:

df.groupby(['sample_id', 'SubSet']).estimate.agg(['mean', 'std']) #.mean()


# In[73]:

# for key, grp in df.groupby('sample_id'):
#     plt.figure()
#     g = sns.boxplot(data=grp, y='SubSet', x='estimate')
#     g.set_title('sample %d' % key)


# What is the "true" cell type?

# In[74]:

test_data['x2']


# In[75]:

#test_df['SubSet'].unique()

true_cell_type = test_data['x2'].transpose()
#true_cell_type[true_cell_type[0] > 0]
true_cell_type


# In[76]:

def savefig(fig, *args, **kwargs):
    """
    Wrap figure.savefig defaulting to tight bounding box.
    From https://github.com/mwaskom/seaborn/blob/dfdd1126626f7ed0fe3737528edecb71346e9eb0/seaborn/axisgrid.py#L1840
    """
    kwargs.setdefault("bbox_inches", "tight")
    fig.savefig(*args, **kwargs)


# In[78]:

for (key, grp), (_, groundtruth) in zip(df.groupby('sample_id'), test_data['x2'].iterrows()):
    f = plt.figure()
    g = sns.boxplot(data=grp, y='SubSet', x='estimate')
    g.set_title('mixture %d' % key)
    
    # add groundtruth points
    gt = groundtruth.copy()
    gt.index = [s.replace('SubSet[', '').replace(']', '') for s in gt.index]
    gt = pd.DataFrame(gt).reset_index()
    gt.columns = ['SubSet', 'estimate']
    #print(gt)
    sns.stripplot(x="estimate", y="SubSet", data=gt,
                  linewidth=0,
                  #jitter=True,
                  #size=3,
                  #color=".3",
                  size=8,
                  color="r",
                  alpha=.8
                 )
    fname = 'plots/oos-mixtures_model6.2_n800MarkerGenes_subset_mixture%d.png' % key
    print(fname)
    savefig(f, fname, dpi=300)
    #break


# In[79]:

# use subplots
f, ax_arr = plt.subplots(5, 2, figsize=(25,25))

for (key, grp), (_, groundtruth), (axid, ax) in zip(df.groupby('sample_id'),
                                                    test_data['x2'].iterrows(),
                                                    np.ndenumerate(ax_arr)):
    #plt.figure()
    g = sns.boxplot(data=grp, y='SubSet', x='estimate', ax=ax)
    g.set_title('mixture %d' % key)
    
    # add groundtruth points
    gt = groundtruth.copy()
    gt.index = [s.replace('SubSet[', '').replace(']', '') for s in gt.index]
    gt = pd.DataFrame(gt).reset_index()
    gt.columns = ['SubSet', 'estimate']
    #print(gt)
    sns.stripplot(x="estimate", y="SubSet", data=gt,
                  ax=ax,
                  linewidth=0,
                  #jitter=True,
                  #size=3,
                  #color=".3",
                  size=8,
                  color="r",
                  alpha=.8
                 )
    #break
#f.show()
savefig(f, 'plots/oos-mixtures_model6.2_n800MarkerGenes_subset_ALL.pdf', dpi=300)


