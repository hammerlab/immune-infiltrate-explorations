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
