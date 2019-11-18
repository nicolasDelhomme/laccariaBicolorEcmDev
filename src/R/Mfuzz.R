#' ---
#' title: "T89 and _Laccaria bicolor_ Mfuzz clustering"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(Mfuzz))

#' * Helpers
source(here("UPSCb-common/src/R/featureSelection.R"))
#source("~/Git/UPSCb-common/src/R/featureSelection.R")

#' * Data
#load(here("analysis/salmon/Potri-all-dds.rda"))
load("/mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor/analysis/salmon/Potri-all-dds.rda")

#' * Function
getCluster <- function(vst,sample.sel,min.std=0.1,clusterN=24){
  
  eset <- ExpressionSet(sapply(split.data.frame(t(vst[,sample.sel]),dds$Time[sample.sel]),colMeans))
  
  f.eset <- filter.std(eset,min.std=min.std)
  
  s.eset <- standardise(f.eset)
  
  m1 <- mestimate(s.eset)
  
  return(list(eset=s.eset,cl=mfuzz(s.eset,c=clusterN,m=m1)))
}

#' # Normalisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#save(vst,file=here("data/analysis/DE/Potri_variance-stabilised-transformed_model-aware_data.rda"))

#' # Clustering
#' ## ECM
ecm <- getCluster(vst,dds$Experiment=="ECM")

#save(ecm$cl,ecm$eset,file="cl.rda")

#' ### Overlapping of the clusters
O <- overlap(ecm$cl)
Ptmp <- overlap.plot(ecm$cl,over=O,thres=0.05)

#' ### Plot clustering results
barplot(table(ecm$cl$cluster),las=2)

mfuzz.plot2(ecm$eset,cl=ecm$cl,mfrow=c(4,4),x1=FALSE)



source("~/Git/UPSCb-common/src/R/gopher.R")

enr.list <- lapply(split(rownames(ecm$eset),ecm$cl$cluster),
    gopher,background=rownames(ecm$eset),
    task="go",url="potri")

#ND: run also for Control (add overlap)
#' ## Control
cont <- getCluster(vst,dds$Experiment=="Cont")

#save(cont$cl,cont$eset,file="cl.rda")

#' ### Overlapping of the clusters
O <- overlap(cont$cl)
Ptmp <- overlap.plot(cont$cl,over=O,thres=0.05)

#' ### Plot clustering results
barplot(table(cont$cl$cluster),las=2)

mfuzz.plot2(cont$eset,cl=cont$cl,mfrow=c(4,4),x1=FALSE)


enr.list.cont <- lapply(split(rownames(cont$eset),cont$cl$cluster),
                   gopher,background=rownames(cont$eset),
                   task="go",url="potri")


#ND: run for ECM - Cont (calculate the ratio ECM - Control)
#' ## ECM - Control
#IS: - ?
#rat$eset <- exprs(ecm$eset) - exprs(cont$eset)

# run comparison of member genes accross the ECM and Cont clusters
# as you said, overlap

#' ## Comparison of member genes accross the ECM and Cont clusters
#' ### ECM's clusters
#' 
#' 
#ecm$cl$membership[ecm$cl$membership > 0.7]

acore.list.ecm <- acore(ecm$eset,cl=ecm$cl,min.acore=0.7)
acore.ecm = list(1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4)

for(i in 1:24){
  acore.ecm[i] <- length(row.names(acore.list.ecm[[i]]))
}

barplot(t(as.matrix(acore.ecm)),las=2,
        main = "Genes belonging to the cluster cores, alpha=0.7, ECM",
        xlab = "Clusters",
        ylab = "Number of genes",
        names.arg = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
        col = "darkred")

#' plot
boxplot(lapply(1:length(ecm$cl$size),function(i){
  ecm$cl$membership[ecm$cl$cluster==i,i]}))

#' ### Control's clusters
acore.list.cont <- acore(cont$eset,cl=cont$cl,min.acore=0.7)

#vector(mode="integer",length=24)
#acore.cont = list(1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4)

#for(i in 1:24){
#  acore.cont[i] <- length(row.names(acore.list.cont[[i]]))
#}

acore.cont <- sapply(acore.list.cont,nrow)

barplot(acore.cont)

barplot(t(as.matrix(acore.cont)),las=2,
        main = "Genes belonging to the cluster cores, alpha=0.7, control",
        xlab = "Clusters",
        ylab = "Number of genes",
        names.arg = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),
        col = "darkred")

# perform the same on the laccaria 
# redo the whole analysis
# take the data from here("data/analysis/DE/Lacbi2....)

load("/mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor/analysis/salmon/Lacbi-all-dds.rda")

# run the full network and use guide genes from T89 and Laccaria
# let's talk when you get there ;-)
# we need a vst with all data, both Potri and LacBi

















