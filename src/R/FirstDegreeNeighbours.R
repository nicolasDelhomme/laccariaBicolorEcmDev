#' ---
#' title: "First Degree Neighbourhood from genes of interest"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' ## Environment
#' library
suppressPackageStartupMessages({
  library(here)
  library(igraph)
  library(SeidRFile)
#  library(pander)
#  library(treemap)
#  library(RColorBrewer)
})

#' # Setup data
#' Read in the backbone network created by seidr at 1%.
sf <- SeidrFileFromPath(here("data/seidr/backbone/backbone-1-percent.sf"))

df <- seidr_iapply(sf, scores = TRUE, ranks = TRUE, 
                  edge_index = TRUE,score_index=which(algorithms(sf) == "irp"),
                   FUN = function(x,nodes){
                     data.frame("From"=nodes[x$node_i], 
                            "To"=nodes[x$node_j])
                     },nodes = nodes(sf))


#' # Import
#' Create an edge-list from the backbone network
graf <- graph.edgelist(as.matrix(df))

#' Check the number of clusters
graf_clust <- clusters(graf)

message(sprintf("There are %s clusters",graf_clust$no))

barplot(table(graf_clust$csize),xlab="number of genes per cluster",ylab="number of clusters")

# Check the proportion of species per cluster
prop <- t(sapply(lapply(
  lapply(split(vertex_attr(graf,"name"),
               graf_clust$membership),substr,1,1),factor,levels=c("P","L")),table))

# Only the first cluster (the largest contains both species)
which(rowSums(prop == 0)<1)

# save the corresponding graph as graphml
dir.create(here("data/analysis/seidr"),showWarnings=FALSE)
write.graph(induced.subgraph(graf, v = get.vertex.attribute(graf,"name")[graf_clust$membership %in% which(rowSums(prop == 0)<1)]),
            file=here("data/analysis/seidr/backbone-1-percent-filtered.graphml"),
            format="graphml")

#' # Gene of interests
#' Load the genes of interests (goi)
gois <- list(
             myb.goi=scan(here("doc/MYB-goi.txt"),what="character"),
             mybr.goi=scan(here("doc/MYB-related-goi.txt"),what="character"),
             wrky.goi=scan(here("doc/WRKY.txt"),what="character"),
             auxin.goi=paste0("Lb",scan(here("doc/auxin-related.txt"),what="character")),
             pectin.goi=paste0("Lb",scan(here("doc/pectin-related.txt"),what="character")))

#' Find the goi present in the network
gois <- lapply(gois,function(goi){goi[goi %in% unlist(df)]})

#' Number of genes present in the network
barplot(sapply(gois,length))

#' Extract the sub-graphs
barons <- lapply(gois,function(goi){
  induced.subgraph(graf, v = goi)
})

#' Remove empty sub-graphs
barons <- barons[sapply(barons,function(baron){
  length(get.vertex.attribute(baron,"name"))
  }) > 0]

#' # First Degree Neighbours (FDN)
#' 
#' ## Extract the X degree neighbours 
#' 
#' (here 1) for every gene, so we get X graphs, one per gene
#' 
#' Which we then reduce and save
subgrafs <- lapply(barons,function(baron){
  Reduce("%u%",make_ego_graph(graf,1,
                 get.vertex.attribute(baron,"name")))
})

dev.null <- sapply(names(subgrafs),function(nam,sg){
  write.graph(graph=sg[[nam]],
              file=here(file.path("data/analysis/seidr",sub("goi","graphml",nam))),
              format="graphml")
},subgrafs)

#' ## Merge
#' Merge the graf and annotate then
#' 
#' Here they are separate but they could be merged and annotated as LD and LL.
fdn.graf <- Reduce("%u%",lapply(names(subgrafs),function(nam,sg){
  set.edge.attribute(sg[[nam]],"dataset",value=sub("\\.goi","",nam))
},subgrafs))

#' Then we create a single dataset attribute and remove the temporary ones
df.network <- paste0(as.integer(!is.na(get.edge.attribute(fdn.network,"dataset_1"))),
                 as.integer(!is.na(get.edge.attribute(fdn.network,"dataset_2"))))

pander(table(df.network))
barplot(table(df.network))

fdn.network <- set.edge.attribute(fdn.network,"dataset", value = df.network)

fdn.network <- remove.edge.attribute(fdn.network,"dataset_1")
fdn.network <- remove.edge.attribute(fdn.network,"dataset_2")
#write_graph(fdn.network, file="/mnt/picea/projects/spruce/meriksson/seidr/FirstDegreeNeighbourhood.graphml",format="graphml")

#' Gain information about the clusters
fdn_clus_LD <- clusters(fdn.graf_LD)
fdn_clus_LL <- clusters(fdn.graf_LL)

#' # Gene Ontology enrichment
#' Get the genes from the networks.
fdn.edges_LD <- as.data.frame(get.edgelist(fdn.graf_LD))
fdn.edges_LL <- as.data.frame(get.edgelist(fdn.graf_LL))

fdn.genes_LD <- unique(unlist(fdn.edges_LD))
fdn.genes_LL <- unique(unlist(fdn.edges_LL))


#' Get the enrichment for the entire networks. Just to see if it makes sense.
source("~/Git/UPSCb/src/R/gopher.R")
fdn_GO_LD <- gopher(fdn.genes_LD,
                    task = list("go", "mapman"),
                    background = bg.unlist_LD, url = "pabies")

fdn_GO_LL <- gopher(fdn.genes_LL,
                    task = list("go", "mapman"),
                    background = bg.unlist_LL, url = "pabies")

#' It's a very global enrichment. For a more in depth enrichment it's better to 
#' look at the individual clusters based on modularity from Gephi.
treemap(fdn_GO_LD$go, 
        index = "name", vSize = "padj", 
        type = "value", vColor = "padj", 
        palette = "PuOr", 
        fontsize.labels=c(15,12), 
        title = "FDN LD - GO",
        inflate.labels = FALSE, 
        lowerbound.cex.labels = 0, 
        bg.labels = "NA", 
        position.legend = "none", 
        file = FALSE)

fdn_GO_LD$mapman$name <- gsub("\\.", " ", fdn_GO_LD$mapman$name)
treemap(fdn_GO_LD$mapman, 
        index = "name", vSize = "padj", 
        type = "value", vColor = "padj", 
        palette = "PuBu", 
        fontsize.labels=c(15,12), 
        title = "FDN LD - mapman", 
        inflate.labels = FALSE, 
        lowerbound.cex.labels = 0, 
        bg.labels = "NA", 
        position.legend = "none", 
        file = FALSE)


#' Once we have obtained the Gephi node/modularity list:
#' We read the exported 'csv' file from Gephi.
#df_gph <- read_delim("MCnodelist.csv", delim = ";", col_names = TRUE, col_types = NULL)

#' Create a background for the gopher enrichment
#merged.bg <- unique(unlist(rbind(bg_LD, bg_LL)))

#' Gopher enrichment of the individual modularity classes 
#lapply(unique(df_gph$modularity_class),
 #      function(mc){enr<-gopher(df_gph$v_names[df_gph$modularity_class == mc],
  #                              task = list("go","mapman"),
   #                             background = bg,
    #                            url = "pabies")
     #  if (length(rownames(enr$go)) > 0) {
      #   go <- enr$go[,c("name","padj")]
       #  file.name <- paste0("results/GO-cluster-", mc,".csv")  
        # write.csv(go, file = file.name, row.names = FALSE)
      # }
      # if (length(rownames(enr$mapman)) > 0) {
      #   mapman <- enr$mapman[,c("name","padj")]
      #   mapman$name <- gsub("\\.", " ", mapman$name)
      #   file.name <- paste0("results/mapman-cluster-", mc,".csv")
      #   write.csv(mapman, file = file.name, row.names = FALSE)  
      # }
       #}) 
#' This function will write 'csv' files with the Gene Ontology (GO) and mapman
#' enrichments is there are any. If there are no enrichments, it will not write a file.
#' The number within the file represents the modularity class and can be traced back
#' to the Gephi network.
