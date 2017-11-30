Now that we are writing paper, looking into provenance of CD8 T cell samples for the "variance is ignored" plot.

# how many samples?

```
 > pyopen column_mappings.pkl

 df=pd.Series(f1).to_frame().reset_index()


 df.columns=['old', 'new']

 In [28]: df.head()
Out[28]:
                              old                            new
0      A_LW_DC48hLPS_HC_U133A.CEL  dendritic | 48hr LPS + HC | 0
1  A_LW_DC48hLPS_U133A_200503.CEL       dendritic | 48hr LPS | 0
2       A_LW_DC6hLPS_HC_U133A.CEL   dendritic | 6hr LPS + HC | 0
3   A_LW_DC6hLPS_U133A_200503.CEL        dendritic | 6hr LPS | 0
4      A_LW_imDC_U133A_200503.CEL       dendritic | Immature | 0


In [24]: df[df['new'].str.contains('CD8')]
Out[24]:
                                                   old                new
62   CD8Tcell-N0-1 [HG-U133A] [IRIS_GSE22886|GSM565...      T cells CD8_1
63   CD8Tcell-N0-2 [HG-U133A] [IRIS_GSE22886|GSM565...      T cells CD8_2
64   CD8Tcell-N0-3 [HG-U133A] [IRIS_GSE22886|GSM565...      T cells CD8_3
65   CD8Tcell-N0-4 [HG-U133A] [IRIS_GSE22886|GSM565...      T cells CD8_4
103                                   GSM155229.CEL.gz    Tcell | CD8 | 0
104                                   GSM155232.CEL.gz    Tcell | CD8 | 1
105                                   GSM155234.CEL.gz    Tcell | CD8 | 2
106                                   GSM155236.CEL.gz    Tcell | CD8 | 3
107                                   GSM155238.CEL.gz    Tcell | CD8 | 4
108                                   GSM565269.CEL.gz  CD8Tcell | N0 | 1
109                                   GSM565270.CEL.gz  CD8Tcell | N0 | 2
110                                   GSM565271.CEL.gz  CD8Tcell | N0 | 3
111                                   GSM565272.CEL.gz  CD8Tcell | N0 | 4



In [27]: print df[df['new'].str.contains('CD8')].old.values
['CD8Tcell-N0-1 [HG-U133A] [IRIS_GSE22886|GSM565269]'
 'CD8Tcell-N0-2 [HG-U133A] [IRIS_GSE22886|GSM565270]'
 'CD8Tcell-N0-3 [HG-U133A] [IRIS_GSE22886|GSM565271]'
 'CD8Tcell-N0-4 [HG-U133A] [IRIS_GSE22886|GSM565272]' 'GSM155229.CEL.gz'
 'GSM155232.CEL.gz' 'GSM155234.CEL.gz' 'GSM155236.CEL.gz'
 'GSM155238.CEL.gz' 'GSM565269.CEL.gz' 'GSM565270.CEL.gz'
 'GSM565271.CEL.gz' 'GSM565272.CEL.gz']

```

matches "Variance within and between cell types" In[59] (in section #1 `T_cell_cols`) 

so 13 samples.


# where are they from?

## master list (`n=9`)


* `CD8Tcell-N0-1 [HG-U133A] [IRIS_GSE22886|GSM565269]`: [Iris](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22886)
* `CD8Tcell-N0-2 [HG-U133A] [IRIS_GSE22886|GSM565270]`: Iris
* `CD8Tcell-N0-3 [HG-U133A] [IRIS_GSE22886|GSM565271]`: Iris
* `CD8Tcell-N0-4 [HG-U133A] [IRIS_GSE22886|GSM565272]`: Iris
* `GSM155229.CEL.gz`: Bindea
* `GSM155232.CEL.gz`: Bindea
* `GSM155234.CEL.gz`: Bindea
* `GSM155236.CEL.gz`: Bindea
* `GSM155238.CEL.gz`: Bindea
* `GSM565269.CEL.gz`: uh oh, a duplicate!
* `GSM565270.CEL.gz`: dupe!
* `GSM565271.CEL.gz`: dupe!
* `GSM565272.CEL.gz`: dupe!

## drop dupes

* `GSM565269.CEL.gz` : `CD8Tcell | N0 | 1`
* `GSM565270.CEL.gz` : `CD8Tcell | N0 | 2`
* `GSM565271.CEL.gz` : `CD8Tcell | N0 | 3`
* `GSM565272.CEL.gz` : `CD8Tcell | N0 | 4`

Made the adjustment in `Variance within and between cell types.ipynb` (Ctrl-F for "drop dupes")

## bindea

in `/curated_data/pure_samples/bindea/cell_type_in_each_sample.tsv`:

```
	type	subtype	filename	filename_low
...
0	Tcell	CD8	GSM155229.CEL.gz	
1	Tcell	CD8	GSM155232.CEL.gz	
2	Tcell	CD8	GSM155234.CEL.gz	
3	Tcell	CD8	GSM155236.CEL.gz	
4	Tcell	CD8	GSM155238.CEL.gz	
...
```


```
> ls /data/microarray/curated_data/pure_samples/bindea/CD8_GSE6740
GSM155229.CEL.gz  GSM155232.CEL.gz  GSM155234.CEL.gz  GSM155236.CEL.gz  GSM155238.CEL.gz
```


Source for those: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6740

> Hyrcza MD, Kovacs C, Loutfy M, Halpenny R et al. Distinct transcriptional profiles in ex vivo CD4+ and CD8+ T cells are established early in human immunodeficiency virus type 1 infection and are characterized by a chronic interferon response as well as extensive transcriptional changes in CD8+ T cells. J Virol 2007 Apr;81(7):3477-86. PMID: 17251300


## iris

Comes from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22886

> Abbas AR, Baldwin D, Ma Y, Ouyang W et al. Immune response in silico (IRIS): immune-specific genes identified from a compendium of microarray expression data. Genes Immun 2005 Jun;6(4):319-31. PMID: 15789058

`curated_data/pure_samples/iris/GSE22886_descriptions_parsed.tsv`:

```
filename	type	subtype	replicate
GSM565269	CD8Tcell	N0	1
GSM565270	CD8Tcell	N0	2
GSM565271	CD8Tcell	N0	3
GSM565272	CD8Tcell	N0	4
```


## cibersort

`curated_data/pure_samples/cibersort/LM22-ref-sample.txt`:

```
# https://github.com/robdmc/pandashells
> cat "curated_data/pure_samples/cibersort/LM22-ref-sample.txt" | p.df -i tsv | p.df 'df.columns'
0
"Relabel"
"A_LW_+BAFF_U133A_250303 [Chtanova_immune|A_LW_+BAFF_U133A_250303]"
"A_LW_+BAFF_U133A_70303 [Chtanova_immune|A_LW_+BAFF_U133A_70303]"
"A_LW_-BAFF_U133A_250303 [Chtanova_immune|A_LW_-BAFF_U133A_250303]"
[...]
"CD8Tcell-N0-1 [HG-U133A] [IRIS_GSE22886|GSM565269]"
"CD8Tcell-N0-2 [HG-U133A] [IRIS_GSE22886|GSM565270]"
"CD8Tcell-N0-3 [HG-U133A] [IRIS_GSE22886|GSM565271]"
"CD8Tcell-N0-4 [HG-U133A] [IRIS_GSE22886|GSM565272]"
[...]
"TH_2 [TREGs_GSE4527|GSM101520]"
"Treg_2 [TREGs_GSE4527|GSM101521]"
```


UH OH we had `CD8Tcell-N0-1 [HG-U133A] [IRIS_GSE22886|GSM565269]` is the same as `GSM565269.CEL.gz`!!! Duplicated!

# how to cite

https://www.ncbi.nlm.nih.gov/geo/info/linking.html

they recommend citing both original paper and GEO accession:

> In March 2008, a search of the GEO Profiles database (Barrett et al., 2006) revealed that Gene X is upregulated in response to compound Y (GEO accession GDSxxx; Smith et al., 2006).

