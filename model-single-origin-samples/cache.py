from stancache.stancache import cached, cached_stan_fit, _read_file
import stancache

stancache.config.set_value(cache_dir = '/modelcache/eliza-immune/cache')
stancache.seed.set_seed()

