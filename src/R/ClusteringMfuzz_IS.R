#' ---
#' title: "T89 and _Laccaria bicolor_ soft clustering with Mfuzz"
#' author: "Nicolas Delhomme && Iryna Shutava"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working directory
setwd("/mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor")
#' ```

#' * Libraries
#suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
# suppressPackageStartupMessages(library(gplots))
# suppressPackageStartupMessages(library(hyperSpec))
# suppressPackageStartupMessages(library(LSD))
# suppressPackageStartupMessages(library(parallel))
# suppressPackageStartupMessages(library(pander))
# suppressPackageStartupMessages(library(plotly))
# suppressPackageStartupMessages(library(RColorBrewer))
# suppressPackageStartupMessages(library(scatterplot3d))
# suppressPackageStartupMessages(library(tidyverse))
# suppressPackageStartupMessages(library(vsn))
library(Mfuzz)
library(cluster)
#library(fpc)

#' * Helper functions
# source("~/Git/UPSCb-common/src/R/plot.multidensity.R")
# source("~/Git/UPSCb-common/src/R/featureSelection.R")

#' * Graphics
# pal <- brewer.pal(8,"Dark2")
# hpal <- colorRampPalette(c("blue","white","red"))(100)
# mar <- par("mar")

load("analysis/salmon/Potri-all-dds.rda")

counts <- counts[,!colnames(counts) %in% c("P11915_104","P11915_127","P11915_139")]

vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#Potri <- new('ExpressionSet', exprs=as.matrix(vst))
sel <- dds$Experiment=="ECM"
eset <- ExpressionSet(sapply(split.data.frame(t(vst[,sel]),dds$Time[sel]),colMeans))





#f.eset <- filter.std(eset,min.std=0.1)
#s.eset <- standardise(f.eset)
#m1 <- mestimate(s.eset)
#cl=mfuzz(s.eset,c=24,m=m1)

Potri <- eset

#' Filtered it. 
#' NB: Works only with "NA", not with "0" value -> in some cases we need to exclude NA genes after standardise
#' 
Potri.r <- filter.NA(Potri, thres=0.25)
#Potri.f <- fill.NA(Potri.r,mode="mean")

Potri.f <- filter.std(Potri.r,min.std=0.1)

Potri.s <- standardise(Potri.f)

#' We need to exclude NA genes.
Potri.s <- filter.NA(Potri.s, thres=0.25)

test <- Potri.s@assayData[["exprs"]]

#' # Soft clustering of gene expression data.
#' ## Initial step - N of the Clusters = 16 
cl <- mfuzz(Potri.s, c=16, m=1.25)

m1 <- mestimate(Potri.s)
m1

cl <- mfuzz(Potri.s, c=16, m=m1)

#ecm <- getCluster(vst,dds$Experiment=="ECM")

#' ### Plot clustering results
#' 
clusplot(exprs(Potri.s), cl$cluster)
points(cl$center,col=1:16,pch=8,cex=1)

#mfuzz.plot2(Potri.s,cl=cl,mfrow=c(4,4),x11=FALSE) # same output as mfuzz.plot
mfuzz.plot2(Potri.s, cl=cl,mfrow=c(4,4),centre=TRUE,x11=FALSE) # lines for cluster centres will be included

#' Extracting list of genes belonging to the cluster cores
list_genes <- acore(Potri.s,cl,min.acore=0.5)

#' ### Overlapping of the clusters
#' 
O <- overlap(cl)
Ptmp <- overlap.plot(cl,over=O,thres=0.05)

#' ## Soft clustering for the N of the clusters = 3
#' 
cl3 <- mfuzz(Potri.s, c=3, m=1.25)

#' ### Plot clustering results
#' 
clusplot(test, cl3$cluster)
points(cl3$center,col=1:16,pch=8,cex=1)

#mfuzz.plot2(Potri.s,cl=cl3,mfrow=c(4,4),x11=FALSE) # same output as mfuzz.plot
mfuzz.plot2(Potri.s, cl=cl3,mfrow=c(4,4),centre=TRUE,x11=FALSE) # lines for cluster centres will be included

#' ### Overlapping of the clusters
#' 
O3 <- overlap(cl3)
Ptmp <- overlap.plot(cl3,over=O3,thres=0.05)




# --------
#' # k-means clustering 
#' 
kl <- kmeans(test, centers = 10, iter.max = 100)
clusplot(test, kl$cluster)
points(kl$center,col=1:10,pch=8,cex=1)

clusplot(test, kl$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

kl2 <- kmeans2(Potri.s,k=16)
clusplot(test, kl2$cluster)
points(kl$center,col=1:16,pch=8,cex=1)


