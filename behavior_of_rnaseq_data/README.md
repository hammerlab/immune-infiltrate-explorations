# How does RNA-seq data behave, in comparison to microarray data?

See Summary notebook for a walk-through of the results.

How this directory is organized:

* `data/`: some datasets downloaded from the RNA-seq processing pipeline I ran
* `data_filenames.tsv`: links data filenames to cell types
* `exploratory/`: some notebooks I was playing with but have not started using
* `plots/`: output plots

Notebooks:

* `0 load.ipynb`: Loads downloaded data and does some processing, e.g. aggregating TPM across genes for each sample. Writes `summary.simple.tsv`. Has some work-in-progress python port from `tximport`.
* `1 EDA.ipynb`: Compares behavior of RNA-seq data to microarray data. Output is in `plots/`.

The data is from [*De novo transcriptome profiling of highly purified human lymphocytes primary cells*](http://www.nature.com/articles/sdata201551), in FASTQ form. I downloaded, trimmed, and ran Kallisto on the data using a Kubernetes cluster. The processing pipeline is in the repo [*hammerlab/infiltrate-rnaseq-pipeline*](https://github.com/hammerlab/infiltrate-rnaseq-pipeline).