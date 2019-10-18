# * Working directory
setwd("/mnt/picea/storage/reference/Populus-trichocarpa/v3.0/")

# * LIbraries
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyverse))
library(seqinr)
library("Biostrings")

tx2gene <- readDNAStringSet("fasta/Ptrichocarpa_v3.0_210_transcript.fa")

tx2gene_tsv = matrix(c(names(tx2gene),names(tx2gene)),nrow = 73013,ncol = 2)
colnames(tx2gene_tsv) <- c("TXID","GENEID")
tx2gene_tsv[,2] <- sub(tx2gene_tsv[,2],pattern = "\\.[0-9]+$",replacement = "")

write.table(tx2gene_tsv, file = "annotation/tx2gene.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
