#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme"
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
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VennDiagram))

#' * Helper files
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds=dds,vst=vst,gene_id){
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    
    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                                 data.frame(value=vst[sel,])),
               aes(x=Time,y=value,col=Experiment,group=Experiment)) +
            geom_point() + geom_smooth() +
            scale_y_continuous(name="VST expression") + 
            ggtitle(label=paste("Expression for: ",gene_id))
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir="analysis/DE",
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds)){
    
    if(length(contrast)==1){
        res <- results(dds,name=contrast)
    } else {
        res <- results(dds,contrast=contrast)
    }
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj)
    
    if(verbose){
        message(sprintf("There are %s genes that are DE",sum(sel)))
    }
            
    if(export){
        if(!dir.exists(default_dir)){
            dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
        }
        write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
        write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
    }
    if(plot){
        heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                  distfun = pearson.dist,
                  hclustfun = function(X){hclust(X,method="ward.D2")},
                  trace="none",col=hpal,labRow = FALSE,
                  labCol=labels[sample_sel]
        )
    }
    return(rownames(res[sel,]))
}

#' # _Laccaria bicolor_
#' * Data
load("analysis/salmon/Lacbi-all-dds.rda")

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Gene of interest
#' * 245379 
#' The gene is medium expressed and show different patterns
#' in ECM and FLM
line_plot(dds,vst,"245379")

#' * 676331
#' The gene is very lowly expressed with some samples having no expression
line_plot(dds,vst,"676331")

#' * 315313
#' Same as 676331
line_plot(dds,vst,"Lacbi1.eu2.Lbscf0069g00950")

#' * 315258
#' The gene has a low expression and there is overrlap between the two
#' experiments, but there is visible expression difference
line_plot(dds,vst,"315258")

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' The model used is:
#' 
#' `Experiment * Time` meaning that the `Experiment` and `Time variable` as 
#' well as their interaction `Experiment:Time` is considered. Because we 
#' cannot assume that these two variables explain all the variance in the data,
#' there is also an `Intercept` for the linear model. This also implies that the 
#' model assumes `ECM` at `3` hours to be the baseline; _i.e._ everything is compared 
#' against it.
resultsNames(dds)

#' ## Results
#' In the following we look at the interaction specific genes; _i.e._ genes that 
#' changes at a given time transition in between experiments
#' ### FLM _vs._ ECM at T3
Lb_3 <- extract_results(dds,vst,"Experiment_FLM_vs_ECM",
                        default_prefix="Labic_FLM-vs-ECM_T3_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==3)

#' ### FLM _vs._ ECM at T7
#' Here we want to conmbine the effect of FLM-ECM at time T3 and the specific
#' FLM:T7 interaction 
Lb_7 <- extract_results(dds,vst,c(0,1,0,0,0,0,1,0,0,0),
                       default_prefix="Labic_FLM-vs-ECM_T7_",
                       labels=paste0(colData(dds)$Experiment,
                                     colData(dds)$Time),
                       sample_sel=colData(dds)$Time==7)

#' ### FLM _vs._ ECM at T14
Lb_14 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,1,0,0),
                   default_prefix="Labic_FLM-vs-ECM_T14_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==14)

#' ### FLM _vs._ ECM at T21
Lb_21 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,1,0),
                        default_prefix="Labic_FLM-vs-ECM_T21_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==21)

#' ### FLM _vs._ ECM at T28
Lb_28 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,0,1),
                        default_prefix="Labic_FLM-vs-ECM_T28_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==28)

#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list(T3=Lb_3,
               T7=Lb_7,
               T14=Lb_14,
               T21=Lb_21,
               T28=Lb_28),
          NULL,
          fill=pal[1:5]))

#' # _Populus tremula_
#' * Data
load("analysis/salmon/Potra-104-127-139-removed-dds.rda")

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Gene of interest
#' * Potra000962g07909
#' The gene is not expressed
line_plot(dds,vst,"Potra000962g07909")

#' * Potra001661g13641
#' The gene is  medium expressed and show time-lagged expression patterns between
#' the experiments
line_plot(dds,vst,"Potra001661g13641")

#' * Potra003711g22520
#' The gene is medium expressed and show initially opposing expression patterns between experiments
line_plot(dds,vst,"Potra003711g22520")

#' * Potra006413g25676
#' The gene his very lowly expressed and shows somewhat opposite expression patterns
line_plot(dds,vst,"Potra006413g25676")

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' The model used is:
#' 
#' `Experiment * Time` meaning that the `Experiment` and `Time variable` as 
#' well as their interaction `Experiment:Time` is considered. Because we 
#' cannot assume that these two variables explain all the variance in the data,
#' there is also an `Intercept` for the linear model. This also implies that the 
#' model assumes `Cont` at `3` hours to be the baseline; _i.e._ everything is compared 
#' against it.
resultsNames(dds)

#' ## Results
#' In the following we look at the interaction specific genes; _i.e._ genes that 
#' changes at a given time transition in between experiments
#' ### ECM _vs._ Cont at T3
Pa_3 <- extract_results(dds,vst,"Experiment_ECM_vs_Cont",
                        default_prefix="Potra_ECM-vs-Cont_T3_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==3)

#' ### ECM _vs._ Cont at T7
Pa_7 <- extract_results(dds,vst,c(0,1,0,0,0,0,1,0,0,0),
                        default_prefix="Potra_ECM-vs-Cont_T7_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==7)

#' ### ECM _vs._ Cont at T14
Pa_14 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,1,0,0),
                         default_prefix="Potra_ECM-vs-Cont_T14_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==14)

#' ### ECM _vs._ Cont at T21
Pa_21 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,1,0),
                         default_prefix="Potra_ECM-vs-Cont_T21_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==21)

#' ### ECM _vs._ Cont at T28
Pa_28 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,0,1),
                         default_prefix="Potra_ECM-vs-Cont_T28_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==28)

#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3,
                            T7=Pa_7,
                            T14=Pa_14,
                            T21=Pa_21,
                            T28=Pa_28),
                       NULL,
                       fill=pal[1:5]))

#' # _Populus trichocarpa_
#' * Data
load("analysis/salmon/Potri-all-dds.rda")

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Gene of interest
goi <- read_lines("~/Git/UPSCb/projects/T89-Laccaria-bicolor/doc/goi.txt")

all(goi %in% rownames(vst))

pdf("candidates-line-plot.pdf",width=12,height=8)
par(mfrow=c(2,3))
lapply(goi,line_plot,dds=dds,vst=vst)
dev.off()

#' * Potri.006G134800
#' The gene is not expressed
line_plot(dds,vst,"Potri.006G134800")

#' * Potri.006G221800
#' The gene is  medium expressed and show time-lagged expression patterns between
#' the experiments
line_plot(dds,vst,"Potri.006G221800")

#' * Potri.002G173900
#' The gene is medium expressed and show initially opposing expression patterns between experiments
line_plot(dds,vst,"Potri.002G173900")

#' * Potri.004G088100
#' The gene his very lowly expressed and shows somewhat opposite expression patterns
line_plot(dds,vst,"Potri.004G088100")

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' The model used is:
#' 
#' `Experiment * Time` meaning that the `Experiment` and `Time variable` as 
#' well as their interaction `Experiment:Time` is considered. Because we 
#' cannot assume that these two variables explain all the variance in the data,
#' there is also an `Intercept` for the linear model. This also implies that the 
#' model assumes `Cont` at `3` hours to be the baseline; _i.e._ everything is compared 
#' against it.
resultsNames(dds)

#' ## Results
#' In the following we look at the interaction specific genes; _i.e._ genes that 
#' changes at a given time transition in between experiments
#' ### ECM _vs._ Cont at T3
Pa_3 <- extract_results(dds,vst,"Experiment_ECM_vs_Cont",
                        default_prefix="Potri_ECM-vs-Cont_T3_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==3)

#' ### ECM _vs._ Cont at T7
Pa_7 <- extract_results(dds,vst,c(0,1,0,0,0,0,1,0,0,0),
                        default_prefix="Potri_ECM-vs-Cont_T7_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==7)

#' ### ECM _vs._ Cont at T14
Pa_14 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,1,0,0),
                         default_prefix="Potri_ECM-vs-Cont_T14_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==14)

#' ### ECM _vs._ Cont at T21
Pa_21 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,1,0),
                         default_prefix="Potri_ECM-vs-Cont_T21_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==21)

#' ### ECM _vs._ Cont at T28
Pa_28 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,0,1),
                         default_prefix="Potri_ECM-vs-Cont_T28_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==28)

#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3,
                            T7=Pa_7,
                            T14=Pa_14,
                            T21=Pa_21,
                            T28=Pa_28),
                       NULL,
                       fill=pal[1:5]))

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


