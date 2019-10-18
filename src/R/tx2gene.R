#' ---
#' title: "Laccaria bicolor transcript to gene"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working directory
setwd("/mnt/picea/storage/reference/Laccaria-bicolor/Lacbi2/gff")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Laccaria-bicolor/Lacbi2/gff")
#' ```

#' * LIbraries
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyverse))

#' # Data
tx2gene <- read_delim("Lacbi2_GeneCatalog_genes_20110203.gff.gz","\t",
           col_names=c("chr","course","type","start","end","score","strand","phase","attribute")) %>% 
    filter(type=="exon") %>% select(value=attribute) %>% 
    mutate(value,value=gsub(value,pattern="name *\"|\"| transcriptId ",
                            replacement="")) %>% 
    separate(1,into=c("GENE","TXID"),sep=";") %>% 
    select("TXID","GENE") %>% distinct()

#' There is only one transcript per gene, or so it seems
ggplot(count(tx2gene,GENE),aes(x=n)) +geom_bar()

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
