import models
import data
import cache
import seed
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import itertools

def execute_models(model_list=['model4', 'model5', 'model6'],
                   by_list=['cell_type', 'SubSet'],
                   sample_n_list=[100, 500],
                   **kwargs):
    
    seed.set_seed()
    exec_params = list(itertools.product(model_list, by_list, sample_n_list))
    results = list()

    for (model, by, sample_n) in exec_params:
        results.append(models.run_model(model_name=model,
                                        sample_n=sample_n,
                                        by=by,
                                        **kwargs))
    return results



