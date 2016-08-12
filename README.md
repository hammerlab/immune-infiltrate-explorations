## Immune infiltrate quantification

We are looking to understand the current state-of-the-art and motivate new methods:

* Understand how existing methods extract representative marker genes or expression profiles from immune cell types, and their drawbacks
* Understand differences between microarry and RNAseq deconvolution
* Understand how existing methods perform deconvolution of an unknown mixture based on reference marker genes or expression signatures, and the drawbacks of these methods
* Show test cases where the existing methods fail
* Collect training/ground truth data from existing papers that we can use
* Prototype new deconvolution models.

See Summary.ipynb for an overview of the analysis. Here's also [the latest project review](https://docs.google.com/presentation/d/1VSIKzc2ygqsOmjYErKSdegzSrqOFhIXRgOGJtcak8VU/edit?usp=sharing).

(Summary.html is a snapshot of Summary.ipynb with better rendering of the table of contents than Github's rendering of Summary.ipynb. It's also saved into gh-pages and hosted at www.hammerlab.org/immune-infiltrate-explorations/Summary.html. The table of contents is made by the TOC jupyter nbextension.)

Directories:

* `curated_data/`: data downloads. Some data is not in Git and is only in a Google bucket -- see README in `curated_data`.
* `data_engineering/`: some transformations of the data
* `evaluate_existing_methods/`: analysis to evaluate existing deconvolution methods
* `hierarchical_models/`: new hierarchical classification models
* `behavior_of_rnaseq_data/`: analysis of how RNA-seq data behaves in comparison to microarray data
* `plots/`: generated figures
* `reference-extraction_writeup/reference_extraction.pdf`: an analysis of the reference extraction -- linked from Summary notebook
