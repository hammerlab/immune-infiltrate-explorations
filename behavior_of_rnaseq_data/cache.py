from stancache.stancache import cached, cached_stan_fit, _read_file
import stancache

stancache.config.set_value(cache_dir = '/mnt/modelcache/immune-infiltrate-explorations')
stancache.seed.set_seed()

