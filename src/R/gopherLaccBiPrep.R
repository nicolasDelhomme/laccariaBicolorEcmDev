#' ---
#' title: "Laccaria bicolor gopher input file generation"
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
  library(parallel)
  library(tidyverse)
  library(here)
})

#' * Out dir
dir.create(here("gopher"),showWarnings = FALSE)

#' # Annotation
#' ## GO
##' GO Black List (Top parents synonyms)
BL <- c("GO:0000004","GO:0007582","GO:0044699","GO:0008150",
        "GO:0005554","GO:0003674","GO:0008372","GO:0005575",
        # and obsolete terms
        "GO:0006803","GO:0006118")

read_tsv(here("fungi/annotation/Lacbi2_GeneCatalog_proteins_20110203_GO.tab.gz"),
               col_types = cols_only(`#proteinId` = col_character(),
                                     goAcc = col_character())) %>% 
  filter(! goAcc %in% BL) %>% 
  rename(ID=`#proteinId`) %>% 
  group_by(ID) %>% summarise(GO=paste(goAcc,collapse="|")) %>% 
  write_tsv(here("gopher/gene_to_go.tsv"),col_names = FALSE)

#' ## KEGG
read_tsv(here("fungi/annotation/Lacbi2_GeneCatalog_proteins_20110203_KEGG.tab.gz"),
         col_types = cols_only(`#proteinId` = col_character(),
                               ecNum = col_character())) %>% 
  rename(ID=`#proteinId`) %>% 
  group_by(ID) %>% summarise(GO=paste(ecNum,collapse="|")) %>% 
  write_tsv(here("gopher/gene_to_kegg.tsv"),col_names = FALSE)

#' ## IPR
read_tsv(here("fungi/annotation/Lacbi2_GeneCatalog_proteins_20110203_IPR.tab.gz"),
         col_types = cols_only(`#proteinId` = col_character(),
                               domainId = col_character(),
                               domainDb = col_character())) %>% 
  filter(domainDb == "HMMPfam") %>%
  rename(ID=`#proteinId`) %>% 
  group_by(ID) %>% summarise(GO=paste(domainId,collapse="|")) %>% 
  write_tsv(here("gopher/gene_to_ipr.tsv"),col_names = FALSE)

#' # Session Info
#' ```{r sessionInfo, echo=FALSE}
#' sessionInfo()
#' ```
