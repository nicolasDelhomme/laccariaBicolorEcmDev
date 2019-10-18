#' ---
#' title: "Tension Wood Biological QA"
#' author: "Nicolas Delhomme & Iryna Shutava"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/jfelten/06_ERF_Project/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/jfelten/06_ERF_Project")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/volcanoPlot.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' Read the sample information
#samples <- read.csv("~/Git/UPSCb/projects/ERF/doc/")

#' # Process
#' Read the HTSeq files in a matrix
res <- mclapply(dir("htseq",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

#' Name them
names(res) <- sub("_sortmerna.*\\.txt","",dir("htseq",pattern="*.txt"))

#' Select
res <- res[grep("_k|_t",names(res))]

#' Create a sample table
samples <- data.frame(Filename=names(res),
                      Sample=sapply(strsplit(names(res),"_"),"[",4))
samples$Tree <- factor(sub("-.*","",samples$Sample))
samples$Section <- factor(sub(".*-","",samples$Sample))

#' Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
res2 <- res
res2[sapply(res,ncol) == 3] <- lapply(res[sapply(res,ncol) == 3],"[",-2)
count.table <- do.call(cbind,lapply(res2,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
dir.create(file.path("analysis","HTSeq"),recursive = TRUE, showWarnings = FALSE)
write.csv(count.table,"analysis/HTSeq/raw-unormalised-data.csv")

#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Convert them into percentages
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' # QA
#' Plot the stats
#' 
#' There are some undersampled libraries
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),names.arg = samples$Sample)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' The average percentage of aligned reads is 56%, which is relatively low,
#' but the sequencing was from very little material and extremely deep
# round(mean(unlist(count.stats["aligned",]/colSums(count.stats))),digits=2)*100
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned",ylim=c(0,1))

#' Check how many genes are never expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative coverage is as expected - around 200X - higher than for
#' standard sequencing project but not so deep due to the large amount of
#' reads aligning in non-genic regions
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by Tissue type (Apical or Suspensor)
#' 
#' The observed distribution is strikingly dissimilar, due to some samples
#' being undersampled
plot.multidensity(log10(count.table),
                  col=pal[as.integer(samples$Tree)],
                  legend.x="topright",
                  legend=levels(samples$Tree),
                  legend.col=pal[1:4],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

plot.multidensity(log10(count.table),
                  col=pal[as.integer(samples$Section)],
                  legend.x="topright",
                  legend=levels(samples$Section),
                  legend.col=pal[1:6],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Remove under-sampled samples
#' all those with less than 1e6 reads
sel <- colSums(count.table) > 1e6
count.table <- count.table[,sel]
samples <- samples[sel,]

#' # Data normalisation 
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' Create the dds object
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(tree=samples$Tree,
                       section=samples$Section),
  design = ~tree)

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation). However, there is again a
#' technical bias, with the aRNA samples being less deep than the original
#' set of samples.
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

#' Validate the VST 
#' 
#' Visualize the corrected mean - sd relationship. It is fairly linear,
#' meaning we can assume homoscedasticity.
#' The slight initial trend / bump is due to genes having few counts in
#' a few subset of the samples and hence having a higher variability. This is
#' expected. 
meanSdPlot(vst[rowSums(count.table)>0,])

#' # QC on the normalised data
#' 
#' ## PCA
#' 
#' First perform a Principal Component Analysis (PCA) of the data
#'to do a quick quality assessment; i.e. replicate should cluster
#' and the first 2-3 dimensions shouldbe explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Plot the PCA first two dimensions
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Tree)],
     pch=c(1:6)[samples$Section],
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",horiz=TRUE,bty=FALSE,
       col=pal[1:nlevels(samples$Tree)],
       pch=16,
       legend = levels(samples$Tree))
legend("bottom",horiz=TRUE,bty=FALSE,
       pch=1:6,
       legend = levels(samples$Section))

#' ## Heatmap
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}

sel <- geneSelect(vst,samples$Sample,exp=3)

dat <- vst[sel,]
heatmap.2(as.matrix(dat),
          labRow=NA,trace="none",
          las=2,col=hpal,labCol = samples$Sample)

#' saturated expression
dat.sat <- dat
plot(density(as.matrix(dat)))
abline(v=c(2,9))
dat.sat[dat < 2 ] <- 2
dat.sat[dat > 9 ] <- 9
heatmap.2(as.matrix(dat.sat),
          labRow=NA,trace="none",
          las=2,col=hpal,
          labCol = samples$Sample)
          
#' standard scores
heatmap.2(as.matrix(dat),
          scale="row",
          labRow=NA,trace="none",
          las=2,col=hpal,
          labCol = samples$Sample)

#' ## Combine the replicates
combined.count.table <- sapply(split.data.frame(t(count.table),samples$Sample),colSums)
write.csv(count.table,"analysis/HTSeq/raw-unormalised_merged-technical-replicates-data.csv")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
