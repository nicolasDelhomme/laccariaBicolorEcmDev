setwd("/mnt/picea/projects/aspseq/jfelten/T89-Laccaria-bicolor")

load("analysis/salmon/Potri-all-dds.rda")

source("~/Git/UPSCb-common/src/R/featureSelection.R")

library(Mfuzz)
library(DESeq2)

vsd <- varianceStabilizingTransformation(dds,blind=FALSE)

vst <- assay(vsd)

vst <- vst - min(vst)

dds.ecm <- dds[,dds$Experiment=="ECM"]
design(dds.ecm) = ~ Time

sel <- featureSelect(vst,dds$Experiment,exp=0.1)

sum(sel)

counts <- vst[sel,]

dat <- sapply(split.data.frame(t(counts),dds$Time),colMeans)
str(dat)

eset <- ExpressionSet(dat)

f.eset <- filter.std(eset,min.std=0.1)

s.eset <- standardise(f.eset)

m1 <- mestimate(s.eset)

m1

cl <- mfuzz(s.eset,c=24,m=m1)

save(cl,s.eset,file="cl.rda")

str(cl)

#mfuzz.plot(s.eset,cl=cl,mfrow=c(4,4),time.labels=seq(0,160,10))

barplot(table(cl$cluster),las=2)

source("~/Git/UPSCb-common/src/R/gopher.R")

enr.list <- lapply(split(rownames(s.eset),cl$cluster),
    gopher,background=rownames(s.eset),
    task="go",url="potri")

# run also for Cont

# run for ECM - Cont

# run comparison of member genes accross the ECM and Cont clusters

# perform the same on the laccaria

# run the full network and use guide genes from T89 and Laccaria



