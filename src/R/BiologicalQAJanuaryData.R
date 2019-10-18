#' ---
#' title: "T89 and _Laccaria bicolor_ Biological QA - January data"
#' author: "Nicolas Delhomme & Iryna Shutava"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working directory
setwd("/mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor/January/Salmon")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor/January/Salmon")
#' ```

#' * Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' * Helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' * Graphics
pal <- brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Metadata
#' Sample information
samples <- read_csv("~/Git/UPSCb/projects/T89-Laccaria-bicolor/doc/Samples.csv") %>% 
  mutate(Time=factor(Time)) %>% 
  mutate(Experiment=factor(Experiment))

#' tx2gene translation
Potra.tx2gene <- read_delim("/mnt/picea/storage/reference/Populus-tremula/v1.1/annotation/tx2gene.tsv",
                            "\t",col_names=c("TXID","GENE"))

Potri.tx2gene <- read_delim("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/annotation/tx2gene.tsv",
                            "\t",col_names=c("TXID","GENE"))

#' # _Laccaria bicolor_ data
#' ## Raw data
lb.filelist <- list.files("Lacbi2", 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)

#' Select the samples containing fungi
stopifnot(all(str_which(basename(lb.filelist),samples$SciLifeID) == 1:nrow(samples)))
names(lb.filelist) <- samples$SciLifeID
lb.filelist <- lb.filelist[samples$Experiment %in% c("ECM","FLM")]

#' Read the expression at the gene level (there is one transcript per gene)
lb.g <- suppressMessages(tximport(files = lb.filelist, 
                                type = "salmon",txOut=TRUE))

counts <- round(lb.g$counts)

#' ## Raw Data QC analysis
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples[match(names(lb.filelist),samples$SciLifeID),]),
       aes(x,y,col=Experiment,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Experiment=samples$Experiment[match(ind,samples$SciLifeID)]) %>% 
  mutate(Time=samples$Time[match(ind,samples$SciLifeID)])

#' Color by Experiment
ggplot(dat,aes(x=values,group=ind,col=Experiment)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' Color by Time
ggplot(dat,aes(x=values,group=ind,col=Time)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(file.path("..","analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file="../analysis/salmon/Lacbi-raw-unormalised-gene-expression_data.csv")

#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
s.sel <- match(colnames(counts),samples$SciLifeID)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples[s.sel,],
  design = ~ Experiment * Time)

save(dds,file="../analysis/salmon/Lacbi-dds.rda")

#' Check the size factors (i.e. the sequencing library size effect)
#' 
#' The sequencing depth is relatively variable (40 to 260 %)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked very well, the data is nearly homoscesdastic
#' 
meanSdPlot(vst[rowSums(counts)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
  
percent <- round(summary(pc)$importance[2,]*100)
  
#' ### 3 first dimensions
#' This looks interesting as the sample separate clearly both by Experiment
#' and Time in the first 2 components.
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Time)[s.sel]],
              pch=c(17,19)[as.integer(samples$Experiment)[s.sel]-1])
legend("topleft",
       fill=pal[1:nlevels(samples$Time)],
       legend=levels(samples$Time))
legend("bottomright",
       pch=c(17,19),
       legend=levels(samples$Experiment)[-1])
par(mar=mar)

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples[s.sel,])

ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

#' ### Heatmap
#' Filter for noise
#' A cutoff at a VST value of 1 leaves about 15000 genes, which is adequate for the QA
conds <- factor(paste(samples$Experiment,samples$Time))[s.sel]
sels <- rangeFeatureSelect(counts=vst,
                   conditions=conds,
                   nrep=3)

#' * Heatmap of "all" genes
#' Taking into account all the genes (above a noise thresholds), the samples cluster
#' according to what we also see in the mapping rate plot, _i.e._ there is a correlation with
#' the amount of sequences in the samples.
heatmap.2(t(scale(t(vst[sels[[2]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

#' *  Heatmap of the 1000 most variable genes
ord <- order(rowSds(vst[sels[[2]],]),decreasing=TRUE) [1:1000]

#' The subset is enriched for higher expression values
ggplot(list(sub=rowMeans(vst[sels[[2]],][ord,]),
            total=rowMeans(vst[sels[[2]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' The clustering remains the same even for the most variable genes
heatmap.2(t(scale(t(vst[sels[[2]],][ord,]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)


#' *  Heatmap of the 1000 least variable genes
ord <- order(rowSds(vst[sels[[2]],])) [1:1000]

#' The subset is enriched for higher expression values, with a strinkingly normal
#' distribution
ggplot(list(sub=rowMeans(vst[sels[[2]],][ord,]),
            total=rowMeans(vst[sels[[2]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' The clustering for the least variable genes shows a separation
#' by experiment and time point
heatmap.2(t(scale(t(vst[sels[[2]],][ord,]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

#' ## Conclusion
#' The quality of the data is good. The PCA shows that the samples cluster by experiment and time, 
#' however, the heatmap shows a clustering that correlates with the mapping rates, _i.e._ the mixed 
#' amount of reads originating from either organism. The final heatmap seem to indicate that this 
#' effect is neglectable albeit confounded.
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' 
#' # T89 ( _Populus tremula_ x _Populus tremuloides_ ) data
#' ## Raw data
pt.filelist <- list.files("Potra01", 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Select the samples containing fungi
stopifnot(all(str_which(basename(pt.filelist),samples$SciLifeID) == 1:nrow(samples)))
names(pt.filelist) <- samples$SciLifeID
pt.filelist <- pt.filelist[samples$Experiment %in% c("Cont","ECM")]

#' Read the expression at the gene level (there is one transcript per gene)
pt.g <- suppressMessages(tximport(files = pt.filelist, 
                                  type = "salmon",
                                  tx2gene=Potra.tx2gene))

counts <- round(pt.g$counts)

#' ## Raw Data QC analysis
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples[match(names(pt.filelist),samples$SciLifeID),]),
       aes(x,y,col=Experiment,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Experiment=samples$Experiment[match(ind,samples$SciLifeID)]) %>% 
  mutate(Time=samples$Time[match(ind,samples$SciLifeID)])

#' Color by Experiment
ggplot(dat,aes(x=values,group=ind,col=Experiment)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' Color by Time
ggplot(dat,aes(x=values,group=ind,col=Time)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(file.path("..","analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file="../analysis/salmon/Potra-raw-unormalised-gene-expression_data.csv")

#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
s.sel <- match(colnames(counts),samples$SciLifeID)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples[s.sel,],
  design = ~ Experiment * Time)

save(dds,file="../analysis/salmon/Potra-dds.rda")

#' Check the size factors (i.e. the sequencing library size effect)
#' 
#' The sequencing depth is highly variable (4 to 900 %)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(counts)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
#' This looks interesting as the sample separate also both by Experiment
#' and Time in the first 2 components, but somewhat not as clearly as for
#' _Laccaria bicolor_
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Time)[s.sel]],
              pch=c(17,19)[as.integer(samples$Experiment)[s.sel]])
legend("topleft",
       fill=pal[1:nlevels(samples$Time)],
       legend=levels(samples$Time))
legend("bottomright",
       pch=c(17,19),
       legend=levels(samples$Experiment)[-3])
par(mar=mar)

#' ### 2D
#' The sample P11915_104 is an outlier
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples[s.sel,])

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

ggplotly(p)

#' ### Heatmap
#' Filter for noise
#' A cutoff at a VST value of 2 leaves about 14000 genes, which is adequate for the QA
conds <- factor(paste(samples$Experiment,samples$Time))[s.sel]
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)

#' * Heatmap of "all" genes
#' Taking into account all the genes (above a noise thresholds), the samples cluster
#' according to what we also see in the mapping rate plot, _i.e._ there is a correlation with
#' the amount of sequences in the samples, but it is not as striking as for _Laccaria bicolor_.
#' There is clustering by time, especially for the earlier time points. The later ECM time 
#' points are closer to the earlier ECM samples, while the later control cluster together, apart from
#'some outliers.
heatmap.2(t(scale(t(vst[sels[[3]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

#' *  Heatmap of the 1000 most variable genes
ord <- order(rowSds(vst[sels[[2]],]),decreasing=TRUE) [1:1000]

#' The subset is enriched for higher expression values
ggplot(list(sub=rowMeans(vst[sels[[2]],][ord,]),
            total=rowMeans(vst[sels[[2]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' The clustering shows a better clustering by time and experiment, even though there are
#' still outliers
heatmap.2(t(scale(t(vst[sels[[2]],][ord,]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)


#' *  Heatmap of the 1000 least variable genes
ord <- order(rowSds(vst[sels[[2]],])) [1:1000]

#' The subset is enriched for higher expression values, with a strinkingly normal
#' distribution
ggplot(list(sub=rowMeans(vst[sels[[2]],][ord,]),
            total=rowMeans(vst[sels[[2]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' The clustering for the least variable genes is very similar to
#' the overall profile
heatmap.2(t(scale(t(vst[sels[[2]],][ord,]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)
#' ## Conclusion
#' The same conclusion as for _Laccaria bicolor_ apply, 
#' despite the difference in the heatmap behaviour that may indicate
#' some bias due to the organisms.
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#----------
#' # _Populus trichocarpa_data
#' ## Raw data
pt.filelist <- list.files("Potri03", 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Select the samples containing fungi
stopifnot(all(str_which(basename(pt.filelist),samples$SciLifeID) == 1:nrow(samples)))
names(pt.filelist) <- samples$SciLifeID
pt.filelist <- pt.filelist[samples$Experiment %in% c("Cont","ECM")]

#' Read the expression at the gene level (there is one transcript per gene)
pt.g <- suppressMessages(tximport(files = pt.filelist, 
                                  type = "salmon",
                                  tx2gene=Potri.tx2gene))

counts <- round(pt.g$counts)

#' ## Raw Data QC analysis
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

ggplot(tibble(x=colnames(counts),y=colSums(counts)) %>% 
         bind_cols(samples[match(names(pt.filelist),samples$SciLifeID),]),
       aes(x,y,col=Experiment,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90),axis.title.x=element_blank())

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Experiment=samples$Experiment[match(ind,samples$SciLifeID)]) %>% 
  mutate(Time=samples$Time[match(ind,samples$SciLifeID)])

#' Color by Experiment
ggplot(dat,aes(x=values,group=ind,col=Experiment)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' Color by Time
ggplot(dat,aes(x=values,group=ind,col=Time)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(file.path("..","analysis","salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file="../analysis/salmon/Potri-raw-unormalised-gene-expression_data.csv")
#' ## Data normalisation 
#' ### Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
s.sel <- match(colnames(counts),samples$SciLifeID)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples[s.sel,],
  design = ~ Experiment * Time)

save(dds,file="../analysis/salmon/Potri-dds.rda")

#' Check the size factors (i.e. the sequencing library size effect)
#' 
#' The sequencing depth is highly variable (4 to 900 %)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(counts)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
#' This looks interesting as the sample separate also both by Experiment
#' and Time in the first 2 components, but somewhat not as clearly as for
#' _Laccaria bicolor_
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Time)[s.sel]],
              pch=c(17,19)[as.integer(samples$Experiment)[s.sel]])
legend("topleft",
       fill=pal[1:nlevels(samples$Time)],
       legend=levels(samples$Time))
legend("bottomright",
       pch=c(17,19),
       legend=levels(samples$Experiment)[-3])
par(mar=mar)

#' ### 2D
#' The sample P11915_104 is an outlier
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples[s.sel,])

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

ggplotly(p)

#' ### Heatmap
#' Filter for noise
#' A cutoff at a VST value of 2 leaves about 14000 genes, which is adequate for the QA
conds <- factor(paste(samples$Experiment,samples$Time))[s.sel]
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)

#' * Heatmap of "all" genes
#' Taking into account all the genes (above a noise thresholds), the samples cluster
#' according to what we also see in the mapping rate plot, _i.e._ there is a correlation with
#' the amount of sequences in the samples, but it is not as striking as for _Laccaria bicolor_.
#' There is clustering by time, especially for the earlier time points. The later ECM time 
#' points are closer to the earlier ECM samples, while the later control cluster together, apart from
#'some outliers.
heatmap.2(t(scale(t(vst[sels[[3]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

#' *  Heatmap of the 1000 most variable genes
ord <- order(rowSds(vst[sels[[2]],]),decreasing=TRUE) [1:1000]

#' The subset is enriched for higher expression values
ggplot(list(sub=rowMeans(vst[sels[[2]],][ord,]),
            total=rowMeans(vst[sels[[2]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' The clustering shows a better clustering by time and experiment, even though there are
#' still outliers
heatmap.2(t(scale(t(vst[sels[[2]],][ord,]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)


#' *  Heatmap of the 1000 least variable genes
ord <- order(rowSds(vst[sels[[2]],])) [1:1000]

#' The subset is enriched for higher expression values, with a strinkingly normal
#' distribution
ggplot(list(sub=rowMeans(vst[sels[[2]],][ord,]),
            total=rowMeans(vst[sels[[2]],])) %>%  melt(),
       aes(x=value,col=L1)) + 
  geom_density() +
  ggtitle("Density of the subset vs. overall") + 
  scale_x_continuous(name=element_text("VST expression")) + 
  theme(legend.title=element_blank())

#' The clustering for the least variable genes is very similar to
#' the overall profile
heatmap.2(t(scale(t(vst[sels[[2]],][ord,]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)
#----------
#' ## Conclusion
#' The same conclusion as for _Laccaria bicolor_ apply, 
#' despite the difference in the heatmap behaviour that may indicate
#' some bias due to the organisms.
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' 
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
