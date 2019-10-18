#' ---
#' title: "T89 and Laccaria bicolor Salmon alignment rates - January data"
#' author: "Nicolas Delhomme"
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
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(tidyverse))

#' # Data
#' ## Metadata
samples <- read_csv("~/Git/UPSCb/projects/T89-Laccaria-bicolor/doc/Samples.csv")

#' ## Lacbi2
lb.files <- list.files("Lacbi2",pattern="*.err",full.names=TRUE)

lb.mapping <- lb.files %>% map_dfr(~read_lines(.x) %>% 
    as_tibble() %>% filter(str_detect(value,pattern="Mapping rate")) %>% 
    mutate(value=gsub(value,pattern=".*= |%",replacement="")) %>% 
    mutate_at(1,as.numeric))

#' ## Potra01
pa.files <- list.files("Potra01",pattern="*.err",full.names=TRUE)

pa.mapping <- pa.files %>% 
    map_dfr(~read_lines(.x) %>% 
                as_tibble() %>% filter(str_detect(value,pattern="Mapping rate")) %>% 
                mutate(value=gsub(value,pattern=".*= |%",replacement="")) %>% 
                mutate_at(1,as.numeric))

#' # Analysis
stopifnot(all(basename(lb.files) == basename(pa.files)))
stopifnot(all(str_which(basename(lb.files),samples$SciLifeID) == 1:nrow(samples)))
dat <- bind_cols(lb.mapping,Potra=pa.mapping) %>% 
    rename_all(~c("Labic","Potra")) %>% 
    bind_cols(samples) %>% mutate(Time=factor(Time))

p <- ggplot(dat,aes(x=Labic,y=Potra,col=Experiment,shape=Time)) + 
    geom_point(size=3) + ggtitle("Mapping rate")

ggplotly(p)

#' # Session Info
#' ```{r session info, eval=FALSE}
#' sessionInfo()
#' ```
