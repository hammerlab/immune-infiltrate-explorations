A good way to look at the results in this directory:

```
mkdir -p mzview
pdftk mixture1.cleaner2.pdf mixture10.cleaner2.pdf mixture2.cleaner2.pdf mixture3.cleaner2.pdf mixture4.cleaner2.pdf mixture5.cleaner2.pdf mixture6.cleaner2.pdf mixture7.cleaner2.pdf mixture8.cleaner2.pdf mixture9.cleaner2.pdf cat output mzview/cleaner2s.pdf
pdftk mixture1.rollup.pdf mixture10.rollup.pdf mixture2.rollup.pdf mixture3.rollup.pdf mixture4.rollup.pdf mixture5.rollup.pdf mixture6.rollup.pdf mixture7.rollup.pdf mixture8.rollup.pdf mixture9.cleaner2.rollup.pdf mixture9.rollup.pdf cat output mzview/rollsups.pdf
pdftk MCSE_sample2-x_dist.pdf Neff_sample2-x_dist.pdf Rhat_sample2-x_dist.pdf StdDev_sample2-x_dist.pdf cat output mzview/diagnostics.pdf
pdftk traceplot_mix10_Treg.pdf traceplot_mix1_memoryB.pdf traceplot_mix1_naiveB.pdf dendrogram.pdf cat output mzview/diagnostics2.pdf
```
