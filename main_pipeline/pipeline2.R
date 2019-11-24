### Pipeline 1.2
# Things accomplished in this pipeline
# Clustering by tree-cut method
# Stacking of clusters for dendrogram visualization
# Idenfication of Hub Genes


######SCRIPT BEGINS HERE

##For one tree, find the dendro of the base...
##Tree Objects
adj_to_tree = function(adj){
  dmat = 1-TOMsimilarity(adj, TOMType="signed")
  geneTree = flashClust(as.dist(dmat), method="average")
  mColorh=NULL
  for (ds in 0:3){
    tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                        minClusterSize = (30-3*ds), cutHeight = 0.99,
                        deepSplit = ds, distM = dmat)
    mColorh=cbind(mColorh,labels2colors(tree$labels));
  }     
  modules = mColorh[,2]
  labels = geneTree$order
  return(list(modules, labels, geneTree))
}

adj_to_order = function(adj){
  dmat = 1-TOMsimilarity(adj, TOMType="signed")
  geneTree = flashClust(as.dist(dmat), method="average")
  return(geneTree$order)
}

##For example, start with 1 - Healthy as base
#Meta function: Compare two adjacency matrices

comparison = function(adjA1, adjA2, str){
  temp = adj_to_tree(adjA1)
  order1 = temp[[2]]
  order2 = adj_to_order(adjB1)
  colors = temp[[1]]
  vec = vector(mode = "character", length = total_genes)
  df = data.frame(label = order1, color = colors)
  for(i in 1:total_genes){
    vec[i] = df[df$label == order2[i],]$color %>% as.character()
  }
  
  mColorh = cbind(colors, vec)
  plotDendroAndColors(temp[[3]], mColorh, c("", ""), main = "", dendroLabels=FALSE);
}

##GSE64913

temp = readRDS("treeObjA3.rds")
order1 = temp[[2]]
order2 = readRDS("treeObjB3.rds")[[3]]$order
colors = temp[[1]]
vec = vector(mode = "character", length = total_genes)
df = data.frame(label = order1, color = colors)
for(i in 1:total_genes){
  vec[i] = df[df$label == order2[i],]$color %>% as.character()
}

mColorh = cbind(colors, vec)
plotDendroAndColors(temp[[3]], mColorh, c("", ""), main = "", dendroLabels=FALSE);

comparison(adjA1, adjB1, "Singhania et al.")
##GSE63142
comparison(adjA2, adjB2, "Wenzel et al.")
##GSE65204
comparison(adjA3, adjB3, "Yang et al.")

triple_comparison = function(adjA1, adjA2, adjA3, str){
  temp = adj_to_tree(adjA1)
  order2 = adj_to_order(adjA2)
  order3 = adj_to_order(adjA3)
  order1 = temp[[2]]
  colors = temp[[1]]
  vec = vector(mode = "character", length = total_genes)
  vec1 = vector(mode = "character", length = total_genes)
  df = data.frame(label = order1, color = colors)
  for(i in 1:total_genes){
    vec[i] = df[df$label == order2[i],]$color %>% as.character()
  }
  for(i in 1:total_genes){
    vec1[i] = df[df$label == order3[i],]$color %>% as.character()
  }
  
  mColorh = cbind(colors, vec, vec1)
  plotDendroAndColors(temp[[3]], mColorh, c("Base", "Comparison 1", "Comparison 2"), main = str,dendroLabels=FALSE);
}
  
triple_comparison(adjA1, adjA2, adjA3, "Clustering Across Healthy")
triple_comparison(adjB1, adjB2, adjB3, "Clustering Across Asthma")


##NEED: Soft Thresholding Data
# Table about clusters obtained for each network
# Hub genes + sub-network analyses: GO, GO, GO, GO v.s. D1, D2, D3 using revigo and such

#For each Module, find the top GO term, (parent GO term) and then do 

n = adj_to_tree(adjA1)[[1]]
w = chooseTopHubInEachModule(t(datExprA1), readRDS("treeObjA1.rds")[[1]], omitColors = F, power = 9)
w = chooseTopHubInEachModule(t(datExprA2), readRDS("treeObjA2.rds")[[1]], omitColors = F, power = 8)
w = chooseTopHubInEachModule(t(datExprA3), readRDS("treeObjA3.rds")[[1]], omitColors = F, power = 7)

w = chooseTopHubInEachModule(t(datExprB1), readRDS("treeObjB1.rds")[[1]], omitColors = F, power = 9)
w = chooseTopHubInEachModule(t(datExprB2), readRDS("treeObjB2.rds")[[1]], omitColors = F, power = 7)
w = chooseTopHubInEachModule(t(datExprB3), readRDS("treeObjB3.rds")[[1]], omitColors = F, power = 7)
