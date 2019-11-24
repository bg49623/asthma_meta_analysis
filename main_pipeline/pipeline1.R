require(ggplot2)
require(dplyr)
require(ggpubr)
require(reshape2)
require(igraph)
require(gtable)
require(DescTools)
require(gridExtra)
require(WGCNA)
source("WGCNA_accessory_functions.R") ##Found in supplementary files
require(flashClust)
require(kableExtra)


process = function(dat){
  # saving gene names
  mapping = as.numeric(dat[, 1])
  # removing gene names from main matrix
  hvals = (dat[, -1])
  # gene names now row names
  rownames(hvals) = mapping
  hvals = as.matrix(hvals)
  hvals = scale(hvals)
  return(hvals)
}
#####LOADING IN STUFF
# A is for healthy, B for asthma
# GSE64913
datExprA1 = process(read.csv("healthy.csv"))
datExprB1 = process(read.csv("disease.csv"))
# GSE63142
datExprA2 = process(read.csv("healthy63.csv"))
datExprB2 = process(read.csv("disease63.csv"))
# GSE65204
datExprA3 = process(read.csv("healthy65.csv"))
datExprB3 = process(read.csv("disease65.csv"))

save(datExprA1, file = "healthy.RData")
save(datExprA2, file = "healthy63.RData")
save(datExprA3, file = "healthy65.RData")
save(datExprB1, file = "disease.RData")
save(datExprB2, file = "disease63.RData")
save(datExprB3, file = "disease65.RData")

load("healthy.RData")
load("healthy63.RData")
load("healthy65.RData")
load("disease.RData")
load("disease63.RData")
load("disease65.RData")

c = goodSamplesGenes(t(datExprA1), verbose = 0)
datExprA1 = datExprA1[c$goodGenes,]
c = goodSamplesGenes(t(datExprA2), verbose = 0)
datExprA2 = datExprA2[c$goodGenes,]
c = goodSamplesGenes(t(datExprA3), verbose = 0)
datExprA3 = datExprA3[c$goodGenes,]
commonProbesA = intersect(intersect (rownames(datExprA1),rownames(datExprA2)), rownames(datExprA3)) 
datExprA1g = datExprA1[commonProbesA,] 
datExprA2g = datExprA2[commonProbesA,] 
datExprA3g = datExprA3[commonProbesA,]

c = NULL
c = goodSamplesGenes(t(datExprB1), verbose = 0)
datExprB1 = datExprB1[c$goodGenes,]
c = goodSamplesGenes(t(datExprB2), verbose = 0)
datExprB2 = datExprB2[c$goodGenes,]
c = goodSamplesGenes(t(datExprA3), verbose = 0)
datExprB3 = datExprB3[c$goodGenes,]
commonProbesA = intersect(intersect (rownames(datExprB1),rownames(datExprB2)), rownames(datExprB3)) 
datExprB1g = datExprB1[commonProbesA,] 
datExprB2g = datExprB2[commonProbesA,] 
datExprB3g = datExprB3[commonProbesA,]

#adjacency matrices
#These take up a lot of memory, apparently; so be conservative when loading them in.

adjA1 = cor(t(datExprA1g))
diag(adjA1) = 0
adjA2 = cor(t(datExprA2g))
diag(adjA2) = 0
adjA3 = cor(t(datExprA3g))
diag(adjA3) = 0
adjB1 = cor(t(datExprB1g))
diag(adjB1) = 0
adjB2 = cor(t(datExprB2g))
diag(adjB2) = 0
adjB3 = cor(t(datExprB3g))
diag(adjB3) = 0

adjA1 = adjacency(t(datExprA1g), type = "signed", power = 9)
diag(adjA1)=0
adjA2 = adjacency(t(datExprA2g), type = "signed", power = 8)
diag(adjA2)=0
adjA3 = adjacency(t(datExprA3g), type = "signed", power = 7)
diag(adjA3)=0
adjB1 = adjacency(t(datExprB1g), type = "signed", power = 9)
diag(adjB1)=0
adjB2 = adjacency(t(datExprB2g), type = "signed", power = 8)
diag(adjB2)=0
adjB3 = adjacency(t(datExprB3g), type = "signed", power = 7)
diag(adjB3)=0

###

#####Clustering Coefficient Graph
cc_graph = function(adj, adj1){
  df = melt(adj)
  df1 = melt(adj1)
    cutoff = 1/10
    df = df[df$value > cutoff,]
    g = graph_from_data_frame(df)
    cc = transitivity(g, type = "local", vids = V(g), isolates = "zero")
    cc = cc[!is.na(cc)]
    df1 = df1[df1$value > cutoff,]
    g = graph_from_data_frame(df1)
    cc1 = transitivity(g, type = "local", vids = V(g), isolates = "zero")
    cc1 = cc1[!is.na(cc1)]
    colors = c(rep("Healthy", length(cc)), rep("Asthma", length(cc1)))
    clustering = c(cc, cc1)
    data.frame(val = clustering, color =  colors) %>%
      ggplot(., aes(x = val, color = color)) + geom_density() + theme_classic() +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      ylab("") + scale_color_brewer(palette="Set1") + theme(legend.position = "none", axis.title.x = element_text(color = "black", size = 12, face = "bold"), panel.border = element_rect(colour = "black", fill=NA, size=1)
) + xlab("Transitivity")
}
pdf("cc1.pdf")
cc_graph(adjA1, adjB1)
dev.off()

pdf("cc2.pdf")
cc_graph(adjA2, adjB2)
dev.off()

pdf("cc3.pdf")
cc_graph(adjA3, adjB3)
dev.off()



ggsave("cc1.pdf", g1, width = 20, height = 3, units = "in", dpi = 300)
ggsave("cc2.pdf", g2, width = 20, height = 3, units = "in", dpi = 300)
ggsave("cc3.pdf", g3, width = 18, height = 3, units = "in", dpi = 300)

print(g1)
print(g2)
print(g3)



ggarrange(g1, g2, g3, nrow = 3)

###Degree Graphs

deg_graph = function(adj){
  df = melt(adj)
  plots = c()
  for(i in 1:9){
    cutoff = i/10
    df = df[df$value > cutoff,]
    g = graph_from_data_frame(df)
    deg = degree(g) %>% as.numeric()
    degs = sort(unique(deg))
    freqs = table(deg) %>% as.numeric()
    plot = data.frame(val = log(degs), freqs = log(freqs)) %>%
      ggplot(., aes(x = val, y = freqs)) + geom_point(shape = 3) + theme_classic() + xlab("Log(Degree)") + 
      ylab("Log(Frequency)") +theme(legend.position = "none", axis.title.x = element_text(color = "black", size = 10, face = "bold"), axis.title.y = element_text(color = "black", size = 10, face = "bold"), panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    plots = c(plots, list(plot))
    P = list(plots=plots, num=9)
  }
  do.call(grid.arrange, c(P$plots, ncol = P$num))
}


g1 = deg_graph(adjA1)
g2 = deg_graph(adjA2)
g3 = deg_graph(adjA3)
g4 = deg_graph(adjB1)
g5 = deg_graph(adjB2)
g6 = deg_graph(adjB3)
ggarrange(g1, g2, g3, g4, g5, g6, nrow = 6)


###Statistics on Network -- Table
# Given an adjacency, find some basic network statistics, then present it in tabular format.
stats = function(adj){
  df = melt(adj)
  df = df[df$value > 0.1,]
  g = graph_from_data_frame(df)
  ##Num. of Edges
  edge = length(E(g))
  ##Average Clustering Coefficient
  ccavg = transitivity(g, type = "local", vids = V(g), isolates = "zero") %>%
    .[!is.na(.)] %>% mean()
  ##Average Degree
  degavg = degree(g) %>% mean()
  ##Average Eigen
  eigavg = eigen_centrality(g)$vector%>% mean()
  ##Maximum k-core
  coremax = coreness(g, mode = "all") %>% max()
  return(c(edge, ccavg, degavg, eigavg, coremax))
}


#save this result, it takes forever to generate

t1 = stats(adjA1)
t2 = stats(adjA2)
t3 = stats(adjA3)
t4 = stats(adjB1)
t5 = stats(adjB2)
t6 = stats(adjB3)

df = as.data.frame(rbind(t1, t2, t3, t4, t5, t6))
colnames(df) = c("# of Edges", "Average Transitivity", "Average Degree", "Average Eigencentrality", "Maximum k-core")
rownames(df) = c("GSE64193-H", "GSE63142-H", "GSE65204-H", "GSE64913-A", "GSE63142-A", "GSE65204-A")
saveRDS(df, "networksummary.rds")

kable(df, "latex", caption = "Group Rows", booktabs = T) %>%
  kable_styling() %>%
  pack_rows("Healthy Networks", 1, 3) %>%
  pack_rows("Asthma Networks", 4, 6)



### Scale Free Plotting

sizeGrWindow(9, 5)
par(mfrow = c(2, 3));
cex1 = 0.9;
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#A1
sft = pickSoftThreshold(t(datExprA1g), powerVector = powers, verbose = 0)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="signed R^2",type="n",
     main = paste("Sighania et al., Healthy"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.81,col="red")

#A2
sft = pickSoftThreshold(t(datExprA2g), powerVector = powers, verbose = 0)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="signed R^2",type="n",
     main = paste("Wenzel et al., Healthy"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.87,col="red")

#A3
sft = pickSoftThreshold(t(datExprA3g), powerVector = powers, verbose = 0)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="signed R^2",type="n",
     main = paste("Yang et al., Healthy"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.93,col="red")

#B1
sft = pickSoftThreshold(t(datExprB1g), powerVector = powers, verbose = 0)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="signed R^2",type="n",
     main = paste("Sighania et al., Asthma"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.72,col="red")

sft = pickSoftThreshold(t(datExprB2g), powerVector = powers, verbose = 0)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="signed R^2",type="n",
     main = paste("Wenzel et al., Asthma"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.83,col="red")

sft = pickSoftThreshold(t(datExprB3g), powerVector = powers, verbose = 0)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="signed R^2",type="n",
     main = paste("Yang et al., Asthma"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")