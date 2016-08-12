Analyses:

* `Compare my LM22 to the paper version.ipynb`: Make matrix with Cibersort and compare
* `Compare basis matrices.ipynb`: Compare and merge LM22 (Cibersort) with Abbas2009 expression basis matrix, and do pairwise Pearson correlations 
* `Compare marker gene lists.ipynb`: See what's in common and what's different between the lists of marker genes for the same cell types.
* `gene_ontology_analysis/`: run GO on the gene lists
* `Evaluating the deconvolution once we have a basis matrix.ipynb`: Test Cibersort deconvolution with synthetic mixtures and find pathological cases. Writes `test_mixtures`.
* `Variance within and between cell types.ipynb`: Look at raw lines before they are made into a basis matrix (e.g. by Cibersort into LM22) to see how much variance there is in expression of all genes and in expression of marker genes, specifically. How much of this information is captured by Cibersort's point estimate?