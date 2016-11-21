
# coding: utf-8

# In[1]:

import data
import models
import cache
import seaborn as sns


# In[2]:

model_name = 'model6.1'
by = 'cell_type'
sample_n = 100


# In[3]:

sample_df = cache.cached(models.prep_sample_df, sample_n=sample_n)
(training_df, test_df) = models.split_sample_df(sample_df=sample_df, test_sample_n=1)


# In[4]:

model_file = models.get_model_file(model_name=model_name)
#print(cache._read_file(model_file))


# In[5]:

stan_data = models.prep_stan_data(sample_df=training_df, test_df=test_df, by=by)


# In[ ]:

model_fit = models.cached_stan_fit(file=model_file, data=stan_data, model_name=model_name)


# In[ ]:



