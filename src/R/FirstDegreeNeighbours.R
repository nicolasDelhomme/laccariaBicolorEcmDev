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
  library(pander)
  library(RColorBrewer)
  library(SeidRFile)
  library(VennDiagram)  
})

#' # Setup 
#' * Graphics
pal <- brewer.pal(8,"Dark2")

#' * Data
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
             wrky.goi=scan(here("doc/WRKY-goi.txt"),what="character"),
             auxin.goi=paste0("Lb",scan(here("doc/auxin-related-goi.txt"),what="character")),
             pectin.goi=paste0("Lb",scan(here("doc/pectin-related-goi.txt"),what="character")))

#' Find the goi present in the network
gois <- lapply(gois,function(goi){goi[goi %in% unlist(df)]})

#' Number of genes present in the network
barplot(sapply(gois,length))

#' Extract the sub-graphs
barons <- lapply(gois,function(goi){
  induced.subgraph(graf, v = goi)
})

#' Remove empty sub-graphs
#barons <- barons[sapply(barons,function(baron){
#  length(get.vertex.attribute(baron,"name"))
#  }) > 0]

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

#' Just a quick look at how many species are found per subgraph
tab <- sapply(lapply(lapply(lapply(subgrafs,get.vertex.attribute,"name"),
                     substr,1,1),factor,levels=c("L","P")),table)
colnames(tab) <- names(subgrafs)
pander(tab)

#' ## Merge
#' Merge the graf and annotate then
#' 
#' Here they are separate but they could be merged and annotated as LD and LL.
fdn.graf <- Reduce("%u%",lapply(names(subgrafs),function(nam,sg){
  attr.nam <- sub("\\.goi","",nam)
  set.edge.attribute(sg[[nam]],attr.nam,value=1)
},subgrafs))

write_graph(fdn.graf, 
            file=here("data/analysis/seidr/firstDegreeNeighbour.graphml"),
            format="graphml")

#' Just a quick look at how many edges are found in common
grid.newpage()
grid.draw(venn.diagram(lapply(edge.attributes(fdn.graf),function(a){which(!is.na(a))}),
                       filename=NULL,fill=pal[1:length(subgrafs)]))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
            