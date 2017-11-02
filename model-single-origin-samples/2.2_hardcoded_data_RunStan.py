
# complements 2.2 Informative priors.ipynb

# to run this interactively while keeping logs:
# {python -i 2.2_RunStan.py} 2>&1 | tee logs/2.2_RunStan.py.consoleout.txt

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

# import informative priors stan data created in the notebook
stan_data = pickle.load(open('2.2_hardcoded-as-data_informative_priors_standata.pkl', 'rb'))
print('managed to load prepared stan data!')



model_name = 'model6.3.1'
by = 'SubSet'

model_file = models.get_model_file(model_name=model_name)
print(model_file)

model_name_for_cache = model_name + '.informative.hardcoded.priors'

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
  pickle.dump(stan_model, open('logs/stanmodel_informative_hardcoded_priors.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)     
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
                    sample_file='logs/informative_hardcoded_priors_sampling_log.txt',
                    #iter=5,
                    #chains=1,
                    #warmup=0
)

end = time()
try:
  elapsed = str(timedelta(seconds=end-start))
  print('{}: Execution completed ({} elapsed)'.format('stan_model.sampling', elapsed))
  pickle.dump(model_fit, open('logs/modelfit_informative_hardcoded_priors.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)     
except:
  print('unexecpted error: {}'.format(sys.exc_info()[0]))

print("Somehow managed to get here without hitting the pystan bug!")


with open('logs/informative_hardcoded_priors_lp.stansummary', 'w') as w:
  w.write(models.filter_stan_summary(model_fit, pars='lp__').to_string())

with open('logs/informative_hardcoded_priors_sample2x.stansummary', 'w') as w:
  w.write(models.filter_stan_summary(model_fit, pars='sample2_x').to_string())
