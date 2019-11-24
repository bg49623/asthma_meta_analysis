### Pipeline 1.4.
# Things done here
# Conversion of edgelist (based on adjacency matrices) into hive plots using the degree framework and HiveR package. 
# Using TFDB and MetaIntegrator results (pipeline 1.0) to create axes categories for hive plots.



require(HiveR)
require(igraph)
require(dplyr)

trr = read.csv("tfdb.csv")
tfs = unique(trr$targetID)
m = read.csv("degenes1.csv", header = FALSE)
library(org.Hs.eg.db)
chars = as.character(m$V1)
omega = mapIds(org.Hs.eg.db, chars, 'ENTREZID', 'SYMBOL')
vec1 = as.numeric(omega)
vec1 = vec1[!is.na(vec1)]
vec = vec1


ugh = readRDS('hive_magenta.rds')
first = ugh[[1]]
g = graph.adjacency(first, weighted = T)
df = get.data.frame(g)
df = df[df$weight > 0.5,]
hnet = graph_from_data_frame(df[,c(1, 2)])
vertices = vertex_attr(hnet)$name %>% as.numeric()
htflist = vertices[vertices %in% tfs] %>% as.numeric()
htarglist = vertices[vertices %in% vec & !(vertices %in% htflist)]
hotherlist = vertices[!(vertices %in% c(htflist,htarglist))]

hhp = edge2HPD(df, axis.cols = c("red", "orange", "yellow"))
hhp$edges = hhp$edges[!(hhp$edges$id1 == hhp$edges$id2), ]
##Extract hive plot nodes dataframe.
hdf = hhp$nodes
##Editing axes
hdf$axis = as.integer(1)
hdf[hdf$lab%in% htarglist, ]$axis = as.integer(2)
hdf[hdf$lab %in% hotherlist, ]$axis = as.integer(3)
########
# CHOOSE AXIS POSITION: Here we use degree
########


dim(hdf)
hdf$radius = degree(hnet, v= V(hnet), mode = "total") %>% jitter()
hdf$color = "white"
hdf$size = 0.5
hhp$nodes = hdf
hhp = manipAxis(hhp, 'norm')
hhp$edges$weight = 0.4

hhpmat = setNames(data.frame(matrix(ncol = 7, nrow = length(hhp$nodes$lab))), c("node.lab", "node.text", "angle", "radius", "offset", "hjust", "vjust"))
library(org.Hs.eg.db)
omega = mapIds(org.Hs.eg.db, hhp$nodes$lab, 'SYMBOL', 'ENTREZID')
hhpmat$node.lab = as.numeric(hhp$nodes$lab)
hhpmat$node.text = as.character(omega)
hhpmat$angle = 270
hhpmat$radius = .5
hhpmat$offset = 0
hhpmat$hjust = 0
hhpmat$vjust = 0
hhpmat$axis = hhp$nodes$axis
hhpmat$pos = hhp$nodes$radius
write.csv(hhpmat, "hhpmat.csv", row.names = FALSE)
plotHive(hhp, axLabs = c("Transcription Factor", "Target Gene", "Other Gene"))


ugh = readRDS('hive_magenta.rds')
first = ugh[[2]]
g = graph.adjacency(first, weighted = T)
df = get.data.frame(g)
df = df[df$weight > 0.7,]
hnet = graph_from_data_frame(df[,c(1, 2)])
vertices = vertex_attr(hnet)$name %>% as.numeric()
htflist = vertices[vertices %in% tfs] %>% as.numeric()
htarglist = vertices[vertices %in% vec & !(vertices %in% htflist)]
hotherlist = vertices[!(vertices %in% c(htflist,htarglist))]

hhp = edge2HPD(df, axis.cols = c("red", "orange", "yellow"))
hhp$edges = hhp$edges[!(hhp$edges$id1 == hhp$edges$id2), ]
##Extract hive plot nodes dataframe.
hdf = hhp$nodes
##Editing axes
hdf$axis = as.integer(1)
hdf[hdf$lab%in% htarglist, ]$axis = as.integer(2)
hdf[hdf$lab %in% hotherlist, ]$axis = as.integer(3)
########
# CHOOSE AXIS POSITION: Here we use degree
########


dim(hdf)
hdf$radius = degree(hnet, v= V(hnet), mode = "total") %>% jitter()
hdf$color = "white"
hdf$size = 0.5
hhp$nodes = hdf
hhp = manipAxis(hhp, 'norm')
hhp$edges$weight = 0.4
plotHive(hhp, axLabs = c("Transcription Factor", "Target Gene", "Other Gene"))


hhpmat = setNames(data.frame(matrix(ncol = 7, nrow = length(hhp$nodes$lab))), c("node.lab", "node.text", "angle", "radius", "offset", "hjust", "vjust"))
library(org.Hs.eg.db)
omega = mapIds(org.Hs.eg.db, hhp$nodes$lab, 'SYMBOL', 'ENTREZID')
hhpmat$node.lab = as.numeric(hhp$nodes$lab)
hhpmat$node.text = as.character(omega)
hhpmat$angle = 270
hhpmat$radius = .5
hhpmat$offset = 0
hhpmat$hjust = 0
hhpmat$vjust = 0
hhpmat$axis = hhp$nodes$axis
hhpmat$pos = hhp$nodes$radius
write.csv(hhpmat, "hhpmat.csv", row.names = FALSE)
plotHive(hhp, axLabs = c("Transcription Factor", "Target Gene", "Other Gene"))





















hiveplot = function(hnet, hubs){

hnet = delete.vertices(simplify(hnet), degree(hnet)==0)
df = as.data.frame(as_edgelist(hnet))
df$value = 1
vertices = vertex_attr(hnet)$name %>% as.numeric()
thehubs = hubs[hubs %in% vertices]

htflist = vertices[vertices %in% tfs] %>% as.numeric()
htarglist = vertices[vertices %in% vec & !(vertices %in% htflist)]
hotherlist = vertices[!(vertices %in% c(htflist,htarglist))]


hhp = edge2HPD(df, axis.cols = c("red", "orange", "yellow"))
hhp$edges = hhp$edges[!(hhp$edges$id1 == hhp$edges$id2), ]
##Extract hive plot nodes dataframe.
hdf = hhp$nodes
##Editing axes
hdf$axis = as.integer(1)
hdf[hdf$lab%in% htarglist, ]$axis = as.integer(2)
hdf[hdf$lab %in% hotherlist, ]$axis = as.integer(3)
########
# CHOOSE AXIS POSITION: Here we use degree
########

dim(hdf)
hdf$radius = degree(hnet, v= V(hnet), mode = "total") %>% jitter()
hdf$color = "white"
hdf$size = 0.5
hdf[hdf$lab %in% thehubs,]$color = "red"
hdf[hdf$lab %in% thehubs,]$size = 1.5
hhp$nodes = hdf
hhp = manipAxis(hhp, 'norm')
hhp$edges$weight = 0.4


hhpmat = setNames(data.frame(matrix(ncol = 7, nrow = length(hhp$nodes$lab))), c("node.lab", "node.text", "angle", "radius", "offset", "hjust", "vjust"))
library(org.Hs.eg.db)
omega = mapIds(org.Hs.eg.db, hhp$nodes$lab, 'SYMBOL', 'ENTREZID')
hhpmat$node.lab = as.numeric(hhp$nodes$lab)
hhpmat$node.text = as.character(omega)
hhpmat$angle = 270
hhpmat$radius = .5
hhpmat$offset = 0
hhpmat$hjust = 0
hhpmat$vjust = 0
hhpmat$axis = hhp$nodes$axis
hhpmat$pos = hhp$nodes$radius
write.csv(hhpmat, "hhpmat.csv", row.names = FALSE)


plotHive(hhp, axLabs = c("Transcription Factor", "Target Gene", "Other Gene"))
}


