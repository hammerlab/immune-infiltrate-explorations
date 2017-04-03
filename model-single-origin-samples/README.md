# How well can we describe RNA-seq data for purified trascript-level data

How this directory is organized:

* `data/`: local copy of datasets downloaded from the RNA-seq processing pipeline, by `data.py`.
* `data_filenames.tsv`: links data filenames to cell types
* `models/`: code for each stan model

See the `SETUP.md` file for description of components necessary to run these on your system.

The data is from [*De novo transcriptome profiling of highly purified human lymphocytes primary cells*](http://www.nature.com/articles/sdata201551), in FASTQ form. I downloaded, trimmed, and ran Kallisto on the data using a Kubernetes cluster. The processing pipeline is in the repo [*hammerlab/infiltrate-rnaseq-pipeline*](https://github.com/hammerlab/infiltrate-rnaseq-pipeline).

Notebooks from exploratory analyses of these data are stored in `../behavior-of-rnaseq`

# Modeling

Modeling in this folder is trying to answer a question:  How well can we describe (purified) RNA-seq transcript-level data using a generative model?


## Process 

There are several variations of models which we have tried to fit to these data. Each model is encoded by a [Stan](http://mc-stan.org) file. In general, I think of each Stan file as a "model", but technically the model (ie likelihood function) is separate from how it's coded.

Model writing & building is iterative, so each model builds on or improves upon the previous model. In some cases, the model improvement is motivated from checking earlier models. In other cases the improvement is to add structure or additional parameters to be estimated. In still other cases the improvement is optimization.


## Stan model files

Model files are stored in the `./models` folder (within this directory). In general, there is at least one workbook for each model, illustrating how well it performs on various (randomly-selected) subsets of our data.


## Sample datasets

Each sample dataset is comprised of data from *all* samples (see EDA notebook) for a subset of genes. In general, the samples are described by the number of genes (n=100, n=500) in the sample. The seed is hard-coded uniformly across all notebooks, so different models estimated for the n=100 sample are estimated on the *same* sample.

Additionally, the sampled data are cached so not only is the seed hard-coded, but the data queried *from the cache* are identical.

## Caching

The caching is handled by the code in the `cache.py` file. All data are pickle files stored in a directory called `.cached_models`, and cached items are named according to the set of keywords (`kwargs`) passed to the function.  Special care was given to things like `dicts` & `pandas.DataFrames` to ensure that the cache name was invariant to the sort order. Thus, you must always used named parameters for cached items.

Because StanFit objects are special (models are compiled first), there is a special caching function to handle these. I've also found they are sensitive to versions of installed software (e.g. of `pystan` and `Cython`), so I have added these to the cache key.

Cached model files can be loaded from the file directly, so long as your currently-installed versions of `pystan` and `Cython` are consistent with the versions at the time the cache was written.

(This probably seems like overkill, but caching these means I can fit the model in one file & then read it in for comparison in another notebook. This keeps things much better organized, believe me!)

