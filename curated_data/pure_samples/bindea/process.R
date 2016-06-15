# per http://bmbolstad.com/misc/ComputeRMAFAQ/ComputeRMAFAQ.html, should process all together
# first, learning to process subdirectory


# http://homer.salk.edu/homer/basicTutorial/affymetrix.html
biocLite('affy')
biocLite('oligo')
biocLite('limma')
library(affy)
list.celfiles()
setwd("C:/Users/maxim/Desktop/code/hammer/infil/curated_data/pure_samples/bindea/garvan")
list.celfiles()
data <- ReadAffy()
eset <- rma(data)
write.exprs(eset,file="log_normalized_expression_mat.txt")



# now processing all.
library(affy)
setwd("C:/Users/maxim/Desktop/code/hammer/infil/curated_data/pure_samples/bindea")
list.celfiles()
data2 <- ReadAffy()
f <- list.files(pattern=".CEL", full.names=T, recursive=T)
f
data2 <- ReadAffy(filenames=f)
eset2 <- rma(data2)
write.exprs(eset2,file="log_normalized_expression_mat_allbindea.txt")
write.exprs(eset2,file="log_normalized_expression_mat_allbindea.csv", sep=',')

# add annotations! like gene name.
biocLite("hgu133a.db")
my_frame <- data.frame(exprs(eset2))
hgu133a()
Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))
hgu113a
hgu113a.db
library("hgu133a.db")
Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))
all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)
head(all)
write.table(all,file="data.ann.txt",sep="\t")
write.table(all,file="data.ann.csv",sep=",")




## ignore below, we are using old array

setwd("C:/Users/maxim/Desktop/code/hammer/infil/curated_data/pure_samples/bindea/garvan")
library(oligo) # note that this overwrites rma method definition above!
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)
library(pd.hg.u133a)
eset_2 <- rma(affyRaw) # this is a different function!! showMethods('rma')
write.exprs(eset_2, file="another_expression_mat.txt")

# some attempts at qc http://master.bioconductor.org/help/course-materials/2009/SeattleApr09/AffyAtoZ/AffymetrixAtoZSlides.pdf
head(sampleNames(eset))
head(sampleNames(data))
library('simpleaffy')
biocLite('simpleaffy')
library('simpleaffy')
saqc <- qc(data)
plot(sac)
plot(saqc)
library("genefilter")
warnings()
dd <- dist2(log2(exprs(data)))
diag(dd) <- 0
dd.row <- as.dendrogram(hclust(as.dist(dd)))
row.ord <- order.dendrogram(dd.row)
library("latticeExtra")
biocLite('latticeExtra')
library("latticeExtra")
legend <- list(top = list(fun = dendrogramGrob,
+ args = list(x = dd.row, side = "top")))
legend <- list(top = list(fun = dendrogramGrob, args = list(x = dd.row, side = "top")))
+ xlab = "", ylab = "", legend = legend)
lp <- levelplot(dd[row.ord, row.ord],scales = list(x = list(rot = 90)),xlab = "", ylab = "", legend = legend)
plot(lp)
library("affyPLM")
biocLite('affyPLM')
library("affyPLM")
dataPLM <- fitPLM(data)
boxplot(dataPLM, main = "NUSE", ylim = c(0.95, 1.22), outline = FALSE, col = "lightblue", las = 3, whisklty = 0, staplelty = 0)
Mbox(dataPLM, main = "RLE", ylim = c(-0.4, 0.4), outline = FALSE, col = "mistyrose", las = 3, whisklty = 0, staplelty = 0)

# another useful read: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf -- e.g. ch6