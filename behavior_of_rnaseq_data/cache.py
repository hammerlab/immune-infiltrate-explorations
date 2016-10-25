import os
import pickle
import pystan
import hashlib
import base64
import logging 
from fnmatch import fnmatch
import ntpath
import seed
from time import time
from datetime import timedelta
import pandas as pd
import re

logger = logging.getLogger(__name__)
_this_dir = os.path.dirname(os.path.abspath(__file__))
cache_dir = os.path.join(_this_dir, '.cached_models')
seed.set_seed()

def _mkdir_if_not_exists(path):
    try:
        os.mkdir(path)
    except:
        pass


def _pickle_dumps_digest(item):
    s = pickle.dumps(item)
    h = _digest(s)
    return h

def _digest(s):
    h = int(hashlib.sha1(s).hexdigest(), 16) % (10 ** 11)
    return h

def _make_digest_dict(k, prefix=''):
    result = dict()
    if len(k) == 0:
        return None
    for (key, item) in sorted(k.items()):
        if isinstance(item, dict):
            item = dict(sorted(item.items()))
            s = _make_digest(item, prefix=key+'-')
            result.update({'{}{}'.format(prefix, key): _digest(s.encode())})
        elif isinstance(item, pd.DataFrame):
            index = tuple(item.index)
            columns = tuple(item.columns)
            values = tuple(tuple(x) for x in item.values)
            s = _pickle_dumps_digest(tuple([index, columns, values]))
            result.update({'{}{}'.format(prefix, key): s})
        else:
            s = _pickle_dumps_digest(item)
            result.update({'{}{}'.format(prefix, key): s})
    return result
    
def _make_digest(k, **kwargs):
    """
    Creates a digest suitable for use within an :class:`phyles.FSCache`
    object from the key object `k`.

    >>> adict = {'a' : {'b':1}, 'f': []}
    >>> make_digest(adict)
    'a2VKynHgDrUIm17r6BQ5QcA5XVmqpNBmiKbZ9kTu0A'
    """
    result = list()
    result_dict = _make_digest_dict(k, **kwargs)
    if result_dict is None:
        return 'default'
    else: 
        for (key, h) in sorted(result_dict.items()):
            result.append('{}_{}'.format(key, h))
        return '.'.join(result)


def _cached_stan_fit(model_name='anon_model', file=None, model_code=None,
                    force=False, cache_dir=cache_dir, cache_only=None, 
                     fit_cachefile=None, **kwargs):
    ''' Cache fit stan model, by storing pickled objects in filesystem
    
    per following warning:
      07: UserWarning: Pickling fit objects is an experimental feature!
        The relevant StanModel instance must be pickled along with this fit object.
        When unpickling the StanModel must be unpickled first.
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
    '''
    ## cached, compiled StanModel
    if file:
        model_code = _read_file(file)
    if model_code:
        model_prefix = '.'.join([model_name, _make_digest(dict(model_code=model_code))])
    if fit_cachefile:
        # if cachefile given, assume cache_only 
        if cache_only is None:
            cache_only = True
        # if necessary, impute cache_dir from filepath
        if fit_cachefile != os.path.basename(fit_cachefile):
            cache_dir, fit_cachefile = os.path.split(os.path.abspath(fit_cachefile))
        # if fit_cachefile given, parse to get fit_model_prefix
        fit_model_prefix = re.sub(string=os.path.basename(fit_cachefile), pattern='(.*_\d+).stanfit.*', repl='\\1')
        if model_code:
            if fit_model_prefix != model_prefix:
                logger.warning('Computed model prefix does not match that used to estimate model. Using prefix matching fit_cachefile')
        model_prefix = fit_model_prefix
    model_cachefile = '.'.join([model_prefix, 'stanmodel', 'pkl'])
    logger.info('Step 1: Get compiled model code, possibly from cache')
    stan_model = cached(func=pystan.StanModel,
                         cache_filename=model_cachefile,
                         model_code=model_code,
                         cache_dir=cache_dir,
                         model_name=model_name,
                         cache_only=cache_only,
                         force=force)
    
    ## cached fitted Stan model, per given args
    if not fit_cachefile:
        fit_cachefile = '.'.join([model_prefix, 'stanfit', _make_digest(dict(**kwargs)), 'pkl'])
    ## either pull fitted model from cache, or fit model
    logger.info('Step 2: Get posterior draws from model, possibly from cache')
    fit = cached(func=stan_model.sampling,
                  cache_filename=fit_cachefile,
                  cache_dir=cache_dir,
                  force=force,
                  cache_only=cache_only,
                  **kwargs)
    logger.info('Fit cachefile: {}'.format(fit_cachefile))
    return fit


def _read_file(filepath):
    with open(filepath, 'r') as myfile:
        data = myfile.read()
    return data


def cached_stan_fit(*args, **kwargs):
    arglist = list(*args)
    if len(arglist)>0:
        raise ValueError('unnamed args not permitted')
    return _cached_stan_fit(**kwargs, seed=seed.seed)


def cached(func, file_prefix='cached', cache_filename=None,
            cache_dir=cache_dir, force=False, cache_only=False,
            compute_hash=True, *args, **kwargs):
    if not cache_filename:
        arglist = list(*args)
        if len(arglist)>0:
            raise ValueError('unnamed args not permitted')
        cache_filename = '.'.join([func.__name__, file_prefix, _make_digest(dict(**kwargs)), 'pkl'])
    logger.info('{}: cache_filename set to {}'.format(func.__name__, cache_filename))
    cache_filepath = os.path.join(cache_dir, cache_filename)
    logger.debug('{}: cache_filepath set to {}'.format(func.__name__, cache_filepath))
    if not force and os.path.exists(cache_filepath):
        try:
            logger.info('{}: Loading result from cache'.format(func.__name__))
            res = pickle.load(open(cache_filepath, 'rb'))
        except:
            print('{}: Error loading from cache'.format(func.__name__))
        else:
            return res
    if cache_only:
        raise ValueError('{}: Cachefile does not exist and cache_only == True. Exiting with failure.'.format(func.__name__))
    logger.info('{}: Starting execution'.format(func.__name__))
    start = time()
    res = func(**kwargs)
    end = time()
    elapsed = str(timedelta(seconds=end-start))
    logger.info('{}: Execution completed ({} elapsed)'.format(func.__name__, elapsed))
    try:
        _mkdir_if_not_exists(cache_dir)
        logger.info('{}: Saving results to cache'.format(func.__name__))
        pickle.dump(res, open(cache_filepath, 'wb'), pickle.HIGHEST_PROTOCOL)     
    except IOError as e:
        logger.warning("{}: I/O error saving to cache ({}): {}".format(func.__name__, e.errno, e.strerror))
    except:
        logger.warning('{}: Unexpected error saving to cache: {}'.format(func.__name__, sys.exc_info()[0]))
    return res

