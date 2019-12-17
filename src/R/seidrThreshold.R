#' ---
#' title: "Seidr Hard Threshold"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' Libs
suppressPackageStartupMessages({
  library(here)
#  library(pander)
  library(reshape2)
#  library(scales)
  library(tidyverse)  
})

#' # Threshold
#' The threshold can be determined by running seidr threshold on the results obtained from 
#' seidr aggregate. The data will be subjected to a threshold range in which each of the genes
#' are assessed on the amount of edges at a specific threshold. The range of the seidr threshold run
#' was between 1 and 0.18
 
#' ## Assessed range
#' The assessed range was between 0.26 and 1.0 because all of the genes were accounted for 
#' at a threshold of 0.26 and only the amount of edges would increase. 
th2 <- read_tsv(file=here("data/seidr/threshold/threshold.txt"),
               col_names=c("Threshold","Edges","Vertices","SFT","ACC")) %>% 
  melt(id = "Threshold")
  
fmt <- function(x){format(x,nsmall = 2, scientific = TRUE)}
ggplot(th2, aes(x = Threshold, y = value, group = variable, col = variable)) +
  geom_line(lwd = 0.5) + facet_wrap(~variable, scales = "free") +
  scale_x_reverse() + theme_bw() +
  theme(text = element_text(size = 10)) +
  scale_y_continuous(labels = fmt)

#' The cutoff can be derived from the obtained graphs above. It is observed that at  
#' a threshold of about 0.4 the SFT plateaus and the ACC start to increaseb exponentially.
#' All edges are accounted for at around 0.35.

#' ## Zoom
#' The graph is zoomed in between a threshold point of 0.25 and 0.45 
#' to devise the cutoff. A final choice is made at 0.3
ggplot(th2, aes(x = Threshold, y = value, group = variable, col = variable)) +
  geom_line(lwd = 0.5) + facet_wrap(~variable, scales = "free") +
  scale_x_reverse() + theme_bw() +
  theme(text = element_text(size = 10)) +
  coord_cartesian(xlim = c(0.25, 0.45)) +
  geom_vline(xintercept=0.3)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
