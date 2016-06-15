Datasets here:

* 
* Bindea with T gamma delta, T follicular helper cells removed (per MSK paper)
*

How to run generators of basis matrices:

* CIBERSORT:
  ```
  java -Xmx3g -Xms3g -jar CIBERSORT.jar -M ExampleMixtures-GEPs.txt -P ~/LM22-ref-sample.txt -c ~/LM22-classes.txt -n 100 -v -q 0.3 -k 999 -m 50 -x 150 -f
True
   ```