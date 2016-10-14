from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pyensembl
from pyensembl import cached_release
import os
from subprocess import call
import tarfile
import pyfscache
import stanity


_this_dir = os.path.dirname(os.path.abspath(__file__))
cache_dir = os.path.join(_this_dir, '.cached_models')
cached = pyfscache.FSCache(cache_dir)

def download(s, data_dir='./data'):
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    ## download & extract data tar file
    data_tar_file = '{}.tar.gz'.format(s)
    call(['gsutil', 'cp', 'gs://mz-hammerlab/output/{}'.format(data_tar_file), data_dir])
    data_tar = tarfile.open(os.path.join(data_dir, data_tar_file))
    data_tar.extractall(path=data_dir)
    data_tar.close()
    ## download & extract log tar file
    log_tar_file = '{}.logs.tar.gz'.format(s)
    call(['gsutil', 'cp', 'gs://mz-hammerlab/output/{}'.format(log_tar_file), data_dir])
    log_tar = tarfile.open(os.path.join(data_dir, data_tar_file))
    log_tar.extractall(path=data_dir)
    log_tar.close()
    print('downloaded', s, 'files to', data_dir)
    os.remove(os.path.join(data_dir, data_tar_file))
    os.remove(os.path.join(data_dir, log_tar_file))


def load_multiple_files(files, ensembl_release=cached_release(79)):
    """
    files is ordered list of samples. they will be assigned sample ids 1-n.
    """
    dfs = []
    for ix, f in enumerate(files):
        filename = 'data/output/%s/abundance.tsv' % f
        if not os.path.exists(filename):
            download(f)
        df = pd.read_csv(filename, sep='\t')
        df['sample_id'] = ix+1
        df['filename'] = f
        df['gene_name'] = df['target_id'].map(lambda t: ensembl_release.gene_name_of_transcript_id(t))
        df['log1p_tpm'] = np.log1p(df['tpm'])
        dfs.append(df)
    return pd.concat(dfs)


def prep_simple_summary(df):
    """
    sum counts, abundance (tpm) across genes
    """
    return df.groupby(['sample_id', 'filename', 'gene_name'])[['est_counts', 'tpm']].sum().reset_index()


## inspired by 
## http://stackoverflow.com/questions/17116814/pandas-how-do-i-split-text-in-a-column-into-multiple-rows
def split_rows_by(df, field, suffix='', by=','):
    s = df[field].str.split(by).apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = '{}{}'.format(field, suffix)
    if suffix == '':
        del df[field]
    return df.join(s)


@cached
def prep_filename_metadata(datafiles=pd.read_csv('data_filenames.tsv', sep='\t')):
    df = datafiles[['filename', 'SubSet', 'Antibody']].copy()
    subsets = split_rows_by(df, 'Antibody')
    subsets['antibody'] = subsets.Antibody.replace(to_replace='.$', value='', inplace=False, regex=True)
    subsets['value'] = subsets.Antibody.replace(to_replace='^.+(.)$', value='\1', inplace=False, regex=True)
    del subsets['Antibody']
    del subsets['SubSet']
    subsets.drop_duplicates(['filename','antibody', 'value'], inplace=True)
    subsets = subsets.pivot('filename', 'antibody', 'value')
    subsets.reset_index(inplace=True)
    subsets = pd.merge(subsets, df.loc[:,['filename','SubSet']], on='filename', how='outer')
    return subsets

@cached
def load_by_cell_type(cell_types, metadata=prep_filename_metadata()):
    file_ids = []
    [file_ids.extend(list(metadata.loc[metadata.SubSet == cell_type, 'filename'].values)) for cell_type in cell_types]
    simple_summary_data = prep_simple_summary(load_multiple_files(file_ids))
    simple_summary_data['log1p_tpm'] = np.log1p(simple_summary_data['tpm'])
    simple_summary_data['log1p_counts'] = np.log1p(simple_summary_data['est_counts'])
    return pd.merge(simple_summary_data, metadata, on='filename', how='left')


def gmean(x):
    return np.expm1(np.mean(np.log1p(x)))


def gstd(x):
    return np.expm1(np.std(np.log1p(x)))


def rescale_geom(x):
    return (x - gmean(x)) / gstd(x) if gstd(x)>0 else 0


@cached
def prep_annotated_data(df):
    df['cell_type'] = df['SubSet'].apply(lambda x: x.split('_')[0])
    df['log1p_tpm_rescaled_type'] = df \
            .groupby(['cell_type','gene_name'])['log1p_tpm'] \
            .transform(rescale_geom)
    df['log1p_tpm_rescaled_subset'] = df \
            .groupby(['SubSet','gene_name'])['log1p_tpm'] \
            .transform(rescale_geom) 
    df['log1p_tpm_rescaled'] = df \
            .groupby(['gene_name'])['log1p_tpm'] \
            .transform(rescale_geom)
    return df

@cached
def cached_stan_fit(*args, **kwargs):
    arglist = list(*args)
    if len(arglist)>0:
        raise ValueError('unnamed args not permitted')
    return stanity.fit(**kwargs)


def print_stan_summary(stan_fit, pars):
    fitsum = stan_fit.summary(pars=pars)
    res = pd.DataFrame(fitsum['summary'], columns=fitsum['summary_colnames'], index=fitsum['summary_rownames'])
    print(res.loc[:,['mean','se_mean','sd','2.5%','50%','97.5%','Rhat']].to_string())

def extract_theta_summary(stan_fit, par='theta', colnames=['B_cell', 'T_cell'], gene_id='new_gene_id'):
    theta = stan_fit.extract(par)[par]
    dflist = list()
    for gene in np.arange(theta.shape[1]):
        part_df = pd.DataFrame(theta[:,gene,:], columns=colnames)
        part_df[gene_id] = gene+1
        part_df.reset_index(inplace=True)
        part_df.rename(columns={'index': 'iter'}, inplace=True)
        dflist.append(part_df)
    return pd.concat(dflist)
    
def prep_theta_summary(stan_fit, sample_df,
                       gene_id='new_gene_id',
                       gene_cat='new_gene_cat',
                       colnames=['B_cell', 'T_cell'],
                       expose_group='T_cell',
                       **kwargs):
    theta_df = extract_theta_summary(stan_fit, gene_id=gene_id, colnames=colnames, **kwargs)
    gene_cat_map = sample_df.drop_duplicates(subset=[gene_id, gene_cat]).loc[:,[gene_id, gene_cat]]
    theta_df = pd.merge(theta_df, gene_cat_map, on=gene_id, how='left')
    theta_ldf = pd.melt(theta_df,
                    value_vars=colnames,
                    value_name='value',
                    id_vars=['iter', gene_id, gene_cat])
    
    ## summarize difference from mean among all items per iteration
    theta_ldf['mean_per_gene'] = theta_ldf.groupby([gene_id, 'iter'])['value'].transform(np.mean)
    theta_ldf['diff_from_mean'] = theta_ldf['value'] - theta_ldf['mean_per_gene']
    ## summarize mean per variable*gene over all iterations
    theta_ldf['mean_value'] = theta_ldf.groupby([gene_id, 'variable'])['value'].transform(np.mean)
    theta_ldf['mean_diff'] = theta_ldf.groupby([gene_id, 'variable'])['diff_from_mean'].transform(np.mean)
    ## expose T-cell value at group level, so this can be sorted on
    if expose_group:
        expose_data = theta_ldf.loc[theta_ldf['variable'] == expose_group, [gene_id, 'mean_value']].drop_duplicates()
        expose_data.rename(columns={'mean_value': '{}_mean'.format(expose_group)}, inplace=True)
        theta_ldf = pd.merge(theta_ldf, expose_data, on=gene_id, how='left')
    return theta_ldf

def prep_sample_df(df, sample_n=None):
    if sample_n:
        sampled_genes = df.drop_duplicates(subset='gene_name').sample(n=sample_n).loc[:,'gene_name']
        sample_df = pd.merge(df, pd.DataFrame(sampled_genes), on='gene_name', how='inner')
    else:
        sample_df = df
    sample_df['new_gene_cat'] = sample_df['gene_name'].astype('category')
    sample_df['new_gene_id'] = sample_df['new_gene_cat'].cat.codes+1
    sample_df['new_sample_cat'] = sample_df['sample_id'].astype('category')
    sample_df['new_sample_id'] = sample_df['new_sample_cat'].cat.codes+1
    return sample_df

def prep_yrep_summary(stan_fit, sample_df, par='y_rep', value_name='pp_est_counts'):
    yrep = stan_fit.extract(par)[par]
    yrep_df = pd.DataFrame(yrep)
    yrep_df.reset_index(inplace=True)
    yrep_df.rename(columns={'index': 'iter'}, inplace=True)
    yrep_ldf = pd.melt(yrep_df, id_vars='iter', value_name=value_name, var_name='index')
    try:
        sample_df.reset_index(inplace=True)
    except:
        pass
    yrep_ldf = pd.merge(yrep_ldf,
                        sample_df,
                        suffixes=['.ppcheck',''],
                        on='index')
    return yrep_ldf


    