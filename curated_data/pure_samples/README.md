Datasets here are labeled as they are cited in the reference extraction writeup (see the md file or bib).

Note that in Bindea dataset, the T gamma delta and T follicular helper cells removed (per Yasin/MSK paper).

The Abbas data is downloaded from the Cibersort website.

Some data is too big for git, and is stored separately.

Data is post-processed in `data_engineering/` directory, where a notebook eventually creates `all_expressions.tsv`.

How to run generators of basis matrices:

* CIBERSORT:
  ```
  java -Xmx3g -Xms3g -jar CIBERSORT.jar -M ExampleMixtures-GEPs.txt -P ~/LM22-ref-sample.txt -c ~/LM22-classes.txt -n 100 -v -q 0.3 -k 999 -m 50 -x 150 -f
True
   ```