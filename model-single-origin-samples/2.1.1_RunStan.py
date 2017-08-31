
# based off 2.1_oos_syntheticmixtures_markergenes.py

# to run this interactively while keeping logs:
# {python -i 2.1.1_RunStan.py} 2>&1 | tee logs/2.1.1_RunStan.py.consoleout.txt

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

# import preprocessed stuff from 2.1.1 notebook
stan_data = pickle.load(open('logs/2.1.1_standata.pkl', 'rb'))
print('managed to load prepared stan data!')



model_name = 'model6.2'
by = 'SubSet'

model_file = models.get_model_file(model_name=model_name)

model_name_for_cache = model_name + '.newmarkergenes'

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
  pickle.dump(stan_model, open('logs/stanmodel_newmarkergenes.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)     
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
                    sample_file='logs/newmarkergenes_sampling_log.txt',
                    #iter=5,
                    #chains=1,
                    #warmup=0
)

end = time()
try:
  elapsed = str(timedelta(seconds=end-start))
  print('{}: Execution completed ({} elapsed)'.format('stan_model.sampling', elapsed))
  pickle.dump(model_fit, open('logs/modelfit_newmarkergenes.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)     
except:
  print('unexecpted error: {}'.format(sys.exc_info()[0]))

print("Somehow managed to get here without hitting the pystan bug!")


with open('logs/markergenes_lp.stansummary', 'w') as w:
  w.write(models.filter_stan_summary(model_fit, pars='lp__').to_string())

with open('logs/markergenes_sample2x.stansummary', 'w') as w:
  w.write(models.filter_stan_summary(model_fit, pars='sample2_x').to_string())