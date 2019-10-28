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

#' * Data
load(here("data/analysis/salmon/Potri-all-dds.rda"))

#' * Function
getCluster <- function(vst,sample.sel,min.std=0.1,clusterN=24){
  
  eset <- ExpressionSet(sapply(split.data.frame(t(vst[,sample.sel]),dds$Time),colMeans))
  
  f.eset <- filter.std(eset,min.std=min.std)
  
  s.eset <- standardise(f.eset)
  
  m1 <- mestimate(s.eset)
  
  return(list(eset=s.eset,cl=mfuzz(s.eset,c=clusterN,m=m1)))
}

#' # Normalisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

save(vst,file=here("data/analysis/DE/Potri_variance-stabilised-transformed_model-aware_data.rda"))

#' # Clustering
#' ## ECM
ecm <- getCluster(vst,dds$Experiment=="ECM")
save(cl,s.eset,file="cl.rda")

#mfuzz.plot(s.eset,cl=cl,mfrow=c(4,4),time.labels=seq(0,160,10))

barplot(table(cl$cluster),las=2)

source("~/Git/UPSCb-common/src/R/gopher.R")

enr.list <- lapply(split(rownames(s.eset),cl$cluster),
    gopher,background=rownames(s.eset),
    task="go",url="potri")

# run also for Control (add overlap)

# run for ECM - Cont (calculate the ratio ECM - Control)
rat$eset <- exprs(ecm$eset) - exprs(cont$eset)

# run comparison of member genes accross the ECM and Cont clusters
# as you said, overlap :-)

# perform the same on the laccaria 
# redo the whole analysis
# take the data from here("data/analysis/DE/Lacbi2....)

# run the full network and use guide genes from T89 and Laccaria
# let's talk when you get there ;-)
# we need a vst with all data, both Potri and LacBi

















