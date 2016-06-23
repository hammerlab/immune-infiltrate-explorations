# http://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/
# http://www.bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsHyperG.pdf
# https://www.bioconductor.org/packages/devel/bioc/manuals/GOstats/man/GOstats.pdf


library(CellMix)
library(BiocInstaller)
biocLite("topGO")
biocLite("GOstats")
biocLite("org.Hs.eg.db")
biocLite("GO.db")
library(org.Hs.eg.db)
library(GO.db)
entrez_object <- org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(entrez_object)
# Convert to a list
entrez_to_go <- as.list(entrez_object[mapped_genes])
#http://www.ncbi.nlm.nih.gov/gene/?term=1
#entrez gene id 1 is A1BG
entrez_to_go[[1]]
go_object <- as.list(org.Hs.egGO2EG)
go_object['GO:0000002']
axon_gene <- go_object['GO:0007411']
length(unlist(axon_gene, use.names=F))
length(unique(unlist(axon_gene, use.names=F)))
axon_gene <- unique(unlist(axon_gene, use.names=F))
head(axon_gene)
library("GOstats")
universe <- mapped_genes
length(axon_gene)
length(universe)
params <- new('GOHyperGParams',
geneIds=axon_gene,
universeGeneIds=universe,
ontology='BP',
pvalueCutoff=0.001,
conditional=F,
testDirection='over',
annotation="org.Hs.eg.db"
)
hgover <- hyperGTest(params)
hgover
result <- summary(hgover)
head(result, 20)


# start
library(CellMix)
cellMarkers()
m <- cellmarkers('IRIS')
details(m)
summary(m)
geneIds(m)
head(geneIds(m))
m$T
geneIds(m)$T
iris_t <- geneIds(m)$T
bindea = read.table('infil/bindea_by_category.csv')
bindea = read.table('infil/bindea_by_category.csv', sep=',')
head(bindea)
bindea = read.table('infil/bindea_by_category.csv', sep=',',rownames=T)
bindea = read.table('infil/bindea_by_category.csv', sep=',',header=T)
head(bindea)
bindea$T.cells
length(bindea$T.cells)
bindea_t = bindea$T.cells
bindea_t = as.character(bindea$T.cells)
bindea_t[lapply(bindea_t, length)>0]
length(bindea_t[lapply(bindea_t, length)>0])
lengths(bindea_t)
bindea_t[bindea_t != ""]
bindea_t = bindea_t[bindea_t != ""]
intersect(iris_t, bindea_t)
set(iris_t) - set(bindea_t)
setdiff(iris_t, bindea_t)
setdiff(bindea_t, iris_t)
s1 = intersect(iris_t, bindea_t)
s2 = setdiff(iris_t, bindea_t)
s3 = setdiff(bindea_t, iris_t)
savehistory("~/testhist.Rhistory")
params <- new('GOHyperGParams',
geneIds=s1,
universeGeneIds=universe,
ontology='BP',
pvalueCutoff=0.001,
conditional=F,
testDirection='over',
annotation="org.Hs.eg.db"
)
universe
params <- new('GOHyperGParams',
geneIds=s1,
universeGeneIds=universe,
ontology='BP',
pvalueCutoff=0.001,
conditional=F,
testDirection='over',
annotation="org.Hs.eg.db"
)


# have to use entrezid
library('annotate')
library('hgu133a.db')
biocLite('hgu133a.db')
library('hgu133a.db')
select(hgu133a.db, s1, "ENTREZID")
select(hgu133a.db, s1, c("SYMBOL", "ENTREZID", "GENENAME"))
select(hgu133a.db, s1, c("SYMBOL", "ENTREZID", "GENENAME"))$ENTREZID
s1e = select(hgu133a.db, s1, c("SYMBOL", "ENTREZID", "GENENAME"))$ENTREZID
s2e = select(hgu133a.db, s2, c("SYMBOL", "ENTREZID", "GENENAME"))$ENTREZID
s3e = select(hgu133a.db, s3, c("SYMBOL", "ENTREZID", "GENENAME"))$ENTREZID
params <- new('GOHyperGParams',
geneIds=s1e,
universeGeneIds=universe,
ontology='BP',
pvalueCutoff=0.001,
conditional=F,
testDirection='over',
annotation="org.Hs.eg.db"
)
hgover <- hyperGTest(params)
hgover
result <- summary(hgover)
head(result, 20)
head(result, 30)
head(result, 50)
head(result, 50)
write.table(result, 'Tcells_intersect_iris_bindea.tsv', sep='\t')
savehistory("~/examplego.Rhistory")
params <- new('GOHyperGParams',
geneIds=s2e,
universeGeneIds=universe,
ontology='BP',
pvalueCutoff=0.001,
conditional=F,
testDirection='over',
annotation="org.Hs.eg.db"
)
hgover <- hyperGTest(params)
result <- summary(hgover)
head(result, 20)
write.table(result, 'Tcells_diff_iris_bindea.tsv', sep='\t')
params <- new('GOHyperGParams',
geneIds=s3e,
universeGeneIds=universe,
ontology='BP',
pvalueCutoff=0.001,
conditional=F,
testDirection='over',
annotation="org.Hs.eg.db"
)
hgover <- hyperGTest(params)
result <- summary(hgover)
head(result, 20)
write.table(result, 'Tcells_diff_bindea_iris.tsv', sep='\t')


# sanity check set behavior
setA<-c("a", "b", "c", "d", "e")
setB<-c("d", "e", "f", "g")
intersect(setA, setB)
setdiff(setA, setB)
setdiff(setB, setA)
