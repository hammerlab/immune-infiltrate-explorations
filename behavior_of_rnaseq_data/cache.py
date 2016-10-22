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

logger = logging.getLogger(__name__)
_this_dir = os.path.dirname(os.path.abspath(__file__))
cache_dir = os.path.join(_this_dir, '.cached_models')
seed.set_seed()

def _mkdir_if_not_exists(path):
    try:
        os.mkdir(path)
    except:
        pass

    
def _make_digest(k):
    """
    Creates a digest suitable for use within an :class:`phyles.FSCache`
    object from the key object `k`.

    >>> adict = {'a' : {'b':1}, 'f': []}
    >>> make_digest(adict)
    'a2VKynHgDrUIm17r6BQ5QcA5XVmqpNBmiKbZ9kTu0A'
    """
    s = pickle.dumps(k)
    h = hashlib.sha256(s).digest()
    b64 = base64.urlsafe_b64encode(h).decode()[:-2]
    return b64.replace('-', '=')


def _cached_stan_fit(model_name='anon_model', file=None, model_code=None,
                    force=False, cache_dir=cache_dir, **kwargs):
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
    model_prefix = '.'.join([model_name, _make_digest(dict(model_code=model_code))])
    model_cachefile = '.'.join([model_prefix, 'stanmodel', 'pkl'])
    logger.info('Step 1: Get compiled model code, possibly from cache')
    stan_model = cached(func=pystan.StanModel,
                         cache_filename=model_cachefile,
                         model_code=model_code,
                         cache_dir=cache_dir,
                         model_name=model_name,
                         force=force)
    
    ## cached fitted Stan model, per given args
    fit_cachefile = '.'.join([model_prefix, 'stanfit', _make_digest(dict(**kwargs)), 'pkl'])
    ## either pull fitted model from cache, or fit model
    logger.info('Step 2: Get posterior draws from model, possibly from cache')
    fit = cached(func=stan_model.sampling,
                  cache_filename=fit_cachefile,
                  cache_dir=cache_dir,
                  force=force,
                  **kwargs)
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
            cache_dir=cache_dir, force=False,
            compute_hash=True, **kwargs):
    if not cache_filename:
        cache_filename = '.'.join([func.__name__, file_prefix, _make_digest(dict(**kwargs)), 'pkl'])
    cache_filepath = os.path.join(cache_dir, cache_filename)
    if not force and os.path.exists(cache_filepath):
        try:
            logger.info('{}: Loading result from cache'.format(func.__name__))
            res = pickle.load(open(cache_filepath, 'rb'))
        except:
            print('{}: Error loading from cache'.format(func.__name__))
        else:
            return res
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

