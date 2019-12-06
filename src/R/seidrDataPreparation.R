#' ---
#' title: "T89 and _Laccaria bicolor_ seidr preparation"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(DESeq2)
  library(plotly)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Metadata
#' Sample information
samples <- read_csv(here("doc/Samples.csv")) %>% 
  mutate(Time=factor(Time)) %>% 
  mutate(Experiment=factor(Experiment))

#' * Laccaria bicolor transcript IDs
tIDs <- scan(here("fungi/fasta/Lacbi2_all_transcripts_20110203-IDs.txt.gz"),
             what="character")

#' * tx2gene
tx2gene <- read_tsv(here("aspen/annotation/tx2gene.tsv"),
                    col_names=TRUE) %>% 
  rbind(data.frame(TXID=tIDs,GENEID=tIDs))


#' # Raw data
filelist <- list.files(here("data/Salmon/Aspen-Laccaria"),
                       recursive = TRUE, 
                       pattern = "quant.sf",
                       full.names = TRUE)

#' Select the samples containing fungi
names(filelist) <- str_match(basename(dirname(filelist)),samples$SciLifeID)

#' Read the expression at the gene level (there is one transcript per gene)
gexp <- suppressMessages(round(tximport(files = filelist,tx2gene=tx2gene,
                                  type = "salmon"))$counts)

#' Sum up by replicates
gexp <- sapply(split.data.frame(t(gexp),colnames(gexp)),colSums)
stopifnot(colnames(gexp)==samples$SciLifeID)

#' ## Data normalisation 
dds <- DESeqDataSetFromMatrix(
  countData = gexp,
  colData = samples,
  design = ~ Experiment * Time)

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked very well, the data is nearly homoscesdastic
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' Definitely not satisfying.
#' 
#' Let's try to process them separately
#' Laccaria
counts <- suppressWarnings(read_csv(here("data/January/analysis/salmon/Lacbi-raw-unormalised-gene-expression_data.csv")) %>% 
                             column_to_rownames("X1") + 
                             read_csv(here("data/analysis/salmon/Lacbi-raw-unormalised-gene-expression_data.csv")) %>% 
                             column_to_rownames("X1"))

s.sel <- match(colnames(counts),samples$SciLifeID)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples[s.sel,],
  design = ~ Experiment * Time)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
lb.vst <- assay(vsd)
rownames(lb.vst)

#' Look OK
meanSdPlot(lb.vst[rowSums(counts)>0,])

#' Aspen
counts <- suppressWarnings(read_csv(here("data/January/analysis/salmon/Potra-raw-unormalised-gene-expression_data.csv")) %>% 
                             column_to_rownames("X1") + 
                             read_csv(here("data/analysis/salmon/Potra-raw-unormalised-gene-expression_data.csv")) %>% 
                             column_to_rownames("X1"))

counts <- counts[,!colnames(counts) %in% c("P11915_104","P11915_127","P11915_139")]
s.sel <- match(colnames(counts),samples$SciLifeID)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples[s.sel,],
  design = ~ Experiment * Time)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pt.vst <- assay(vsd)
#' Look OK
meanSdPlot(pt.vst[rowSums(counts)>0,])

plot(density(lb.vst),xlim=c(0,15),col=2,lwd=2,lty=2)

lines(density(pt.vst),col=3,lty=3)

pt.vst <- pt.vst + (min(lb.vst) - min(pt.vst))

lines(density(pt.vst),col=4,lty=3)

str(lb.vst)

vst <- matrix(min(pt.vst),ncol=nrow(samples),nrow=sum(nrow(lb.vst),nrow(pt.vst)),
              dimnames=list(c(rownames(lb.vst),rownames(pt.vst)),samples$SciLifeID))

vst[rownames(lb.vst),colnames(lb.vst)] <- lb.vst
vst[rownames(pt.vst),colnames(pt.vst)] <- pt.vst

samples <- samples[match(colnames(vst),samples$SciLifeID),]

vst <- vst - min(vst)

lines(density(vst),col=1,lty=2)

pc <- prcomp(t(as.data.frame(vst)))

percent <- round(summary(pc)$importance[2,]*100)

pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples)

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Experiment,text=SciLifeID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(name=element_text(paste("PC1 (",percent[1],"%)",sep=""))) +
  scale_y_continuous(name=element_text(paste("PC2 (",percent[2],"%)",sep="")))

ggplotly(p)

vst <- vst[,!colnames(vst) %in% c("P11915_104","P11915_127","P11915_139")]

#' # Filter
sels <- rangeFeatureSelect(as.matrix(vst),
                          conditions=as.factor(paste(
                            samples$Experiment,
                            samples$Time)),nrep=3)
cutoff <- 1

rownames(vst) <- sub("\\|.*","",sub("jgi\\|Lacbi2\\|","Lb",rownames(vst)))


#' # Export
dir.create(here("data/seidr"),showWarnings=FALSE)

#' * gene by column, without names matrix
write.table(t(vst[sels[[cutoff+1]],]),
            file=here("data/seidr/headless.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(sub("\\.1$","",rownames(vst)[sels[[cutoff+1]]])),
            file=here("data/seidr/genes.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)


#' ## Conclusion
#' The quality of the data is good. The PCA shows that the samples cluster by experiment and time, 
#' however, the heatmap shows a clustering that correlates with the mapping rates, _i.e._ the mixed 
#' amount of reads originating from either organism. The final heatmap seem to indicate that this 
#' effect is neglectable albeit confounded.
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' 
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
