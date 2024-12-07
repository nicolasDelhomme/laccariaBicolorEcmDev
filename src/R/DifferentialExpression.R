#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme & Iryna Shutava"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' Load the genes of interests (goi)
gois <- list(
    myb.goi=scan(here("doc/MYB-goi.txt"),what="character"),
    mybr.goi=scan(here("doc/MYB-related-goi.txt"),what="character"),
    wrky.goi=scan(here("doc/WRKY-goi.txt"),what="character"),
    auxin.goi=scan(here("doc/auxin-related-goi.txt"),what="character"),
    pectin.goi=scan(here("doc/pectin-related-goi.txt"),what="character"))

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
                              export=TRUE,default_dir=here("data/analysis/DE"),
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
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' # _Laccaria bicolor_
#' * Data
load(here("data/analysis/salmon/Lacbi-all-dds.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Genes from the GOI
rownames(vst) <- substr(rownames(vst),12,17)

#' * Genes from the GOI
lapply(gois, function(goi){
        pdf(here("data/analysis/DE",paste0(goi,"-Lacbi-candidates-line-plot.pdf")),width=12,height=8)
        par(mfrow=c(2,3))
        lapply(goi,function(x){
        if ((x %in% rownames(vst))){
            line_plot(dds,vst,x)
        }
    })
    dev.off()
})

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
                        default_prefix="Lacbi_FLM-vs-ECM_T3_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==3)

#' ### FLM _vs._ ECM at T7
#' Here we want to conmbine the effect of FLM-ECM at time T3 and the specific
#' FLM:T7 interaction 
Lb_7 <- extract_results(dds,vst,c(0,1,0,0,0,0,1,0,0,0),
                       default_prefix="Lacbi_FLM-vs-ECM_T7_",
                       labels=paste0(colData(dds)$Experiment,
                                     colData(dds)$Time),
                       sample_sel=colData(dds)$Time==7)

#' ### FLM _vs._ ECM at T14
Lb_14 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,1,0,0),
                   default_prefix="Lacbi_FLM-vs-ECM_T14_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==14)

#' ### FLM _vs._ ECM at T21
Lb_21 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,1,0),
                        default_prefix="Lacbi_FLM-vs-ECM_T21_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==21)

#' ### FLM _vs._ ECM at T28
Lb_28 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,0,1),
                        default_prefix="Lacbi_FLM-vs-ECM_T28_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==28)

#' ### Venn Diagram
#' #### All DE genes
grid.newpage()
grid.draw(venn.diagram(list(`3d`=Lb_3$all,
                            `7d`=Lb_7$all,
                            `14d`=Lb_14$all,
                            `21d`=Lb_21$all,
                            `28d`=Lb_28$all),
                       NULL,
                       fill=pal[1:5]))

pdf(file=here("data/analysis/DE/Labic_FLM-vs-ECM_VennDiagram-DE-all.pdf"),
    width=10,height=10)
grid.newpage()
grid.draw(venn.diagram(list(`3d`=Lb_3$all,
                            `7d`=Lb_7$all,
                            `14d`=Lb_14$all,
                            `21d`=Lb_21$all,
                            `28d`=Lb_28$all),
                       NULL,
                       fill=pal[1:5]))
dev.off()

#' Export the common and unique sets
write(file=here("data/analysis/DE/Labic_FLM-vs-ECM_common.txt"),
    Reduce("intersect",list(`3d`=Lb_3$all,
                        `7d`=Lb_7$all,
                        `14d`=Lb_14$all,
                        `21d`=Lb_21$all,
                        `28d`=Lb_28$all)))

write(file=here("data/analysis/DE/Labic_FLM-vs-ECM_3d.txt"),
      setdiff(Lb_3$all,Reduce("union",list(Lb_7$all,
                          Lb_14$all,
                          Lb_21$all,
                          Lb_28$all))))

write(file=here("data/analysis/DE/Labic_FLM-vs-ECM_7d.txt"),
      setdiff(Lb_7$all,Reduce("union",list(Lb_3$all,
                                           Lb_14$all,
                                           Lb_21$all,
                                           Lb_28$all))))

write(file=here("data/analysis/DE/Labic_FLM-vs-ECM_14d.txt"),
      setdiff(Lb_14$all,Reduce("union",list(Lb_7$all,
                                           Lb_3$all,
                                           Lb_21$all,
                                           Lb_28$all))))

write(file=here("data/analysis/DE/Labic_FLM-vs-ECM_21d.txt"),
      setdiff(Lb_21$all,Reduce("union",list(Lb_7$all,
                                           Lb_14$all,
                                           Lb_3$all,
                                           Lb_28$all))))

write(file=here("data/analysis/DE/Labic_FLM-vs-ECM_28d.txt"),
      setdiff(Lb_28$all,Reduce("union",list(Lb_7$all,
                                           Lb_14$all,
                                           Lb_21$all,
                                           Lb_3$all))))

#' #### UP DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Lb_3$up,
                            T7=Lb_7$up,
                            T14=Lb_14$up,
                            T21=Lb_21$up,
                            T28=Lb_28$up),
                       NULL,
                       fill=pal[1:5]))

pdf(file=here("data/analysis/DE/Labic_FLM-vs-ECM_VennDiagram-DE-up.pdf"),
    width=10,height=10)
grid.newpage()
grid.draw(venn.diagram(list(`3d`=Lb_3$up,
                            `7d`=Lb_7$up,
                            `14d`=Lb_14$up,
                            `21d`=Lb_21$up,
                            `28d`=Lb_28$up),
                       NULL,
                       fill=pal[1:5]))
dev.off()

#' #### DOWN DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Lb_3$dn,
                            T7=Lb_7$dn,
                            T14=Lb_14$dn,
                            T21=Lb_21$dn,
                            T28=Lb_28$dn),
                       NULL,
                       fill=pal[1:5]))

pdf(file=here("data/analysis/DE/Labic_FLM-vs-ECM_VennDiagram-DE-dn.pdf"),
    width=10,height=10)
grid.newpage()
grid.draw(venn.diagram(list(`3d`=Lb_3$dn,
                            `7d`=Lb_7$dn,
                            `14d`=Lb_14$dn,
                            `21d`=Lb_21$dn,
                            `28d`=Lb_28$dn),
                       NULL,
                       fill=pal[1:5]))
dev.off()

res.list <- list(Lb_3=list(all=substr(Lb_3$all,12,17),
                           up=substr(Lb_3$up,12,17),
                           dn=substr(Lb_3$dn,12,17)),
                 Lb_7=list(all=substr(Lb_7$all,12,17),
                           up=substr(Lb_7$up,12,17),
                           dn=substr(Lb_7$dn,12,17)),
                 Lb_14=list(all=substr(Lb_14$all,12,17),
                            up=substr(Lb_14$up,12,17),
                            dn=substr(Lb_14$dn,12,17)),
                 Lb_21=list(all=substr(Lb_21$all,12,17),
                            up=substr(Lb_21$up,12,17),
                            dn=substr(Lb_21$dn,12,17)),
                 Lb_28=list(all=substr(Lb_28$all,12,17),
                            up=substr(Lb_28$up,12,17),
                            dn=substr(Lb_28$dn,12,17)))

#' ### Gene Ontology enrichment
#' ```{r go, echo=FALSE,eval=FALSE}
#' Once you have obtained a list of candidate genes, you most probably want
#' to annotate them.
#' 
#' In the following example, we first identify the background; _i.e._ the
#' population of expressed genes. We select the genes expressed in a least
#' 2 replicate of one condition at a cutoff of `exp`.
#' 
#' Next we run the enrichment, in the example against `athaliana` using 
#' the gofer3 REST API (interfaced through the gopher.R script loaded at the
#' beginning of this fil).
#' 
#' Finally we export the go enrichment as a complete table and as a table consisting
#' of only the `id` and `padj` columns. The latter can be used as input for _e.g._
#' REVIGO.
#' ```
background <- rownames(vst)[featureSelect(vst,dds$Experiment,exp=0.1)]

enr.list <- lapply(res.list,function(r){
    lapply(r,gopher,background=background,task="go",url="lacbi2")
})

dev.null <- lapply(names(enr.list),function(n){
    r <- enr.list[[n]]
    if (! is.null(r$all$go)){
        write_delim(r$all$go,path=file.path(file.path(here("data/analysis/DE",
                                                       paste0(n,"-all-DE-genes_GO-enrichment.txt")))))
        write_delim(r$all$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                        paste0(n,"-all-DE-genes_GO-enrichment_for-REVIGO.txt")))))
    }
    if (! is.null(r$up$go)){
        write_csv(r$up$go,path=file.path(file.path(here("data/analysis/DE",
                                                    paste0(n,"-up-DE-genes_GO-enrichment.txt")))))
        write_delim(r$up$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                       paste0(n,"-up-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
    }
    if (! is.null(r$dn$go)){
        write_csv(r$dn$go,path=file.path(file.path(here("data/analysis/DE",
                                                    paste0(n,"-down-DE-genes_GO-enrichment.txt")))))
        write_delim(r$dn$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                       paste0(n,"-down-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
    }
})


#' # _Populus tremula_
#' * Data
load(here("data/analysis/salmon/Potra-104-127-139-removed-dds.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Gene of interest
#' * Potra000962g07909
#' The gene is not expressed
line_plot(dds,vst,"Potra000962g07909")

#' * Genes from the GOI
lapply(gois, function(goi){
    pdf(here("data/analysis/DE",paste0(goi,"-Potra-candidates-line-plot.pdf")),width=12,height=8)
    par(mfrow=c(2,3))
    lapply(goi,function(x){
        if ((x %in% rownames(vst))){
            line_plot(dds,vst,x)
        }
    })
    dev.off()
})

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

#' #### Extra validation
r.sel <- rownames(vst) %in% Pa_7$up
c.sel <- dds$Time == "7"

dat <- t(scale(t(vst[r.sel,]))) %>% 
  as.data.frame() %>% 
  rownames_to_column("GeneID") %>% 
  pivot_longer(starts_with("P11915")) %>% 
  left_join(samples %>% dplyr::select(SciLifeID,Experiment,Time),by = c("name"="SciLifeID"))

ggplot(dat,
       aes(x=Time,y=value,group=c(GeneID))) +
  geom_smooth(se=FALSE) +
  scale_y_continuous(name="VST expression") +
  facet_wrap(~Experiment)

ggplot(dat,
       aes(x=parse_integer(as.character(Time)),
           y=value,group=Time,fill=Time)) +
  #geom_boxplot() +
  geom_violin(trim=FALSE,draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_continuous(name="standard score") +
  scale_x_continuous(paste0("days (n=",length(unique(dat$GeneID)),")")) +
  facet_wrap(~Experiment)

r.sel <- rownames(vst) %in% Pa_7$dn
dat <- t(scale(t(vst[r.sel,]))) %>% 
  as.data.frame() %>% 
  rownames_to_column("GeneID") %>% 
  pivot_longer(starts_with("P11915")) %>% 
  left_join(samples %>% dplyr::select(SciLifeID,Experiment,Time),by = c("name"="SciLifeID"))

ggplot(dat,
       aes(x=Time,y=value,group=c(GeneID))) +
  geom_smooth(se=FALSE) +
  scale_y_continuous(name="VST expression") +
  facet_wrap(~Experiment)

ggplot(dat,
       aes(x=parse_integer(as.character(Time)),
           y=value,group=Time,fill=Time)) +
  #geom_boxplot() +
  geom_violin(trim=FALSE,draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_continuous(name="standard score") +
  scale_x_continuous(paste0("days (n=",length(unique(dat$GeneID)),")")) +
  facet_wrap(~Experiment)

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
#' #### All DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3$all,
                            T7=Pa_7$all,
                            T14=Pa_14$all,
                            T21=Pa_21$all,
                            T28=Pa_28$all),
                       NULL,
                       fill=pal[1:5]))

pdf(file=here("data/analysis/DE/Potra_ECM-vs-Cont_VennDiagram-DE-all.pdf"),
    width=10,height=10)
grid.newpage()
grid.draw(venn.diagram(list(`3d`=Pa_3$all,
                            `7d`=Pa_7$all,
                            `14d`=Pa_14$all,
                            `21d`=Pa_21$all,
                            `28d`=Pa_28$all),
                       NULL,
                       fill=pal[1:5]))
dev.off()

#' Export the common and unique sets
write(file=here("data/analysis/DE/Potra_ECM-vs-Cont_common.txt"),
      Reduce("intersect",list(`3d`=Pa_3$all,
                          `7d`=Pa_7$all,
                          `14d`=Pa_14$all,
                          `21d`=Pa_21$all,
                          `28d`=Pa_28$all)))

write(file=here("data/analysis/DE/Potra_ECM-vs-Cont_3d.txt"),
      setdiff(Pa_3$all,Reduce("union",list(Pa_7$all,
                                           Pa_14$all,
                                           Pa_21$all,
                                           Pa_28$all))))

write(file=here("data/analysis/DE/Potra_ECM-vs-Cont_7d.txt"),
      setdiff(Pa_7$all,Reduce("union",list(Pa_3$all,
                                           Pa_14$all,
                                           Pa_21$all,
                                           Pa_28$all))))

write(file=here("data/analysis/DE/Potra_ECM-vs-Cont_14d.txt"),
      setdiff(Pa_14$all,Reduce("union",list(Pa_7$all,
                                            Pa_3$all,
                                            Pa_21$all,
                                            Pa_28$all))))

write(file=here("data/analysis/DE/Potra_ECM-vs-Cont_21d.txt"),
      setdiff(Pa_21$all,Reduce("union",list(Pa_7$all,
                                            Pa_14$all,
                                            Pa_3$all,
                                            Pa_28$all))))

write(file=here("data/analysis/DE/Potra_ECM-vs-Cont_28d.txt"),
      setdiff(Pa_28$all,Reduce("union",list(Pa_7$all,
                                            Pa_14$all,
                                            Pa_21$all,
                                            Pa_3$all))))



#' #### UP DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3$up,
                            T7=Pa_7$up,
                            T14=Pa_14$up,
                            T21=Pa_21$up,
                            T28=Pa_28$up),
                       NULL,
                       fill=pal[1:5]))

pdf(file=here("data/analysis/DE/Potra_ECM-vs-Cont_VennDiagram-DE-up.pdf"),
    width=10,height=10)
grid.newpage()
grid.draw(venn.diagram(list(`3d`=Pa_3$up,
                            `7d`=Pa_7$up,
                            `14d`=Pa_14$up,
                            `21d`=Pa_21$up,
                            `28d`=Pa_28$up),
                       NULL,
                       fill=pal[1:5]))
dev.off()

#' #### DOWN DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3$dn,
                            T7=Pa_7$dn,
                            T14=Pa_14$dn,
                            T21=Pa_21$dn,
                            T28=Pa_28$dn),
                       NULL,
                       fill=pal[1:5]))

pdf(file=here("data/analysis/DE/Potra_ECM-vs-Cont_VennDiagram-DE-dn.pdf"),
    width=10,height=10)
grid.newpage()
grid.draw(venn.diagram(list(`3d`=Pa_3$dn,
                            `7d`=Pa_7$dn,
                            `14d`=Pa_14$dn,
                            `21d`=Pa_21$dn,
                            `28d`=Pa_28$dn),
                       NULL,
                       fill=pal[1:5]))
dev.off()

res.list <- list(Pa_3=Pa_3,
                 Pa_7=Pa_7,
                 Pa_14=Pa_14,
                 Pa_21=Pa_21,
                 Pa_28=Pa_28)

#' ### Gene Ontology enrichment
#' ```{r go_2, echo=FALSE,eval=FALSE}
#' Once you have obtained a list of candidate genes, you most probably want
#' to annotate them.
#' 
#' In the following example, we first identify the background; _i.e._ the
#' population of expressed genes. We select the genes expressed in a least
#' 2 replicate of one condition at a cutoff of `exp`.
#' 
#' Next we run the enrichment, in the example against `athaliana` using 
#' the gofer3 REST API (interfaced through the gopher.R script loaded at the
#' beginning of this fil).
#' 
#' Finally we export the go enrichment as a complete table and as a table consisting
#' of only the `id` and `padj` columns. The latter can be used as input for _e.g._
#' REVIGO.
#' ```
background <- rownames(vst)[featureSelect(vst,dds$Experiment,exp=0.1)]

enr.list <- lapply(res.list,function(r){
    lapply(r,gopher,background=background,task="go",url="potra")
})

dev.null <- lapply(names(enr.list),function(n){
    r <- enr.list[[n]]
    if (! is.null(r$all$go)){
        write_delim(r$all$go,path=file.path(file.path(here("data/analysis/DE",
                                                           paste0(n,"-all-DE-genes_GO-enrichment.txt")))))
        write_delim(r$all$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                            paste0(n,"-all-DE-genes_GO-enrichment_for-REVIGO.txt")))))
    }
    if (! is.null(r$up$go)){
        write_csv(r$up$go,path=file.path(file.path(here("data/analysis/DE",
                                                        paste0(n,"-up-DE-genes_GO-enrichment.txt")))))
        write_delim(r$up$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                           paste0(n,"-up-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
    }
    if (! is.null(r$dn$go)){
        write_csv(r$dn$go,path=file.path(file.path(here("data/analysis/DE",
                                                        paste0(n,"-down-DE-genes_GO-enrichment.txt")))))
        write_delim(r$dn$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/DE",
                                                                           paste0(n,"-down-DE-genes_GO-enrichment_for-REVIGO.txt")))))    
    }
})

#' # _Populus trichocarpa_
#' * Data
load(here("data/analysis/salmon/Potri-all-dds.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

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
#' #### All DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3$all,
                            T7=Pa_7$all,
                            T14=Pa_14$all,
                            T21=Pa_21$all,
                            T28=Pa_28$all),
                       NULL,
                       fill=pal[1:5]))

#' #### UP DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3$up,
                            T7=Pa_7$up,
                            T14=Pa_14$up,
                            T21=Pa_21$up,
                            T28=Pa_28$up),
                       NULL,
                       fill=pal[1:5]))
#' #### DOWN DE genes
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3$dn,
                            T7=Pa_7$dn,
                            T14=Pa_14$dn,
                            T21=Pa_21$dn,
                            T28=Pa_28$dn),
                       NULL,
                       fill=pal[1:5]))

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


