# make hierarchical clustering plots (heatmap w/ dendrogram) of LM22 and Abbas data
# this is how plots/lm22.pdf was made.
# (plots/lm22.png was made from plots/lm22.pdf by imagemagick's convert command in Summary.ipynb)
# run the below, then manually save off the graphics


setwd('C:\\Users\\maxim\\Desktop\\code\\hammer\\infil')
biocLite('ALL')
library('ALL')
lm22=read.table('LM22.txt', sep='\t', header=T, row.names='Gene.symbol')
head(lm22)
heatmap(as.matrix(lm22))

abbas=read.table('GSE11103_matrix_pure.txt', sep='\t', header=T, row.names='X.Sample_title')
head(abbas)
heatmap(as.matrix(abbas))


heatmap
dist
hclust
