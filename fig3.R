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

## Singhania et al.
adj1 = cor(t(datExprA1g))
adj2 = cor(t(datExprB1g))
adj = abs(adj2 - adj1)

scores = c()
cutoffs = seq(0.5, 2, by = 0.1)
for(i in 5:20){
      cutoff = i/10
      adj[adj < cutoff] = NA
      df = melt(adj, na.rm = TRUE)
      df = df[df$Var1 %in% vec,]
      g = graph_from_data_frame(df)
      denn1 = V(g)$name
      scores = c(scores, length(denn1))
}

data.frame(val = cutoffs, freqs = scores) %>% ggplot(., aes(x = val, y = freqs)) + geom_point(shape = 20, size = 5) + theme_bw(base_size = 24) + xlab("Difference Cutoff") + ylab("DE Genes and Nearest Neighbors")+ theme(legend.position = "none", axis.title.x = element_text(color = "black", size = 16, face = "bold"), axis.title.y = element_text(color = "black", size = 16, face = "bold"), panel.border = element_rect(colour = "black", fill=NA, size=1)) + geom_vline(xintercept = 0.9, size = 1.5, linetype = "dashed", color = "red") + main("title");


adj[adj < 1.3] = NA #lower if needed
df = melt(adj, na.rm = TRUE)
df = df[df$Var1 %in% vec,]
g1 = graph_from_data_frame(df)
denn1 = V(g1)$name

## Wenzel et al.
adj1 = cor(t(datExprA2g))
adj2 = cor(t(datExprB2g))
adj = abs(adj2 - adj1)

scores = c()
cutoffs = seq(0.5, 2, by = 0.1)
for(i in 5:20){
      cutoff = i/10
      adj[adj < cutoff] = NA
      df = melt(adj, na.rm = TRUE)
      df = df[df$Var1 %in% vec,]
      g = graph_from_data_frame(df)
      denn1 = V(g)$name
      scores = c(scores, length(denn1))
}

data.frame(val = cutoffs, freqs = scores) %>%
      ggplot(., aes(x = val, y = freqs)) + geom_point(shape = 3, size = 2) + theme_classic() + xlab("Difference Cutoff") + 
      ylab("DE Genes and Nearest Neighbors")+ theme(legend.position = "none", axis.title.x = element_text(color = "black", size = 10, face = "bold"), axis.title.y = element_text(color = "black", size = 10, face = "bold"), panel.border = element_rect(colour = "black", fill=NA, size=1));

adj[adj < 0.9] = NA #lower if needed
df = melt(adj, na.rm = TRUE)
df = df[df$Var1 %in% vec,]
g2 = graph_from_data_frame(df)
denn2 = V(g2)$name


## Yang et al. 

adj1 = cor(t(datExprA3g))
adj2 = cor(t(datExprB3g))
adj = abs(adj2 - adj1)


scores = c()
cutoffs = seq(0.5, 2, by = 0.1)
for(i in 5:20){
      cutoff = i/10
      adj[adj < cutoff] = NA
      df = melt(adj, na.rm = TRUE)
      df = df[df$Var1 %in% vec,]
      g = graph_from_data_frame(df)
      denn1 = V(g)$name
      scores = c(scores, length(denn1))
}

data.frame(val = cutoffs, freqs = scores) %>%
      ggplot(., aes(x = val, y = freqs)) + geom_point(shape = 3, size = 2) + theme_classic() + xlab("Difference Cutoff") + 
      ylab("DE Genes and Nearest Neighbors")+ theme(legend.position = "none", axis.title.x = element_text(color = "black", size = 10, face = "bold"), axis.title.y = element_text(color = "black", size = 10, face = "bold"), panel.border = element_rect(colour = "black", fill=NA, size=1));

adj[adj < 0.9] = NA #lower if needed
df = melt(adj, na.rm = TRUE)
df = df[df$Var1 %in% vec,]
g3 = graph_from_data_frame(df)
denn3 = V(g3)$name

denn = denn1 %u% denn2 %u% denn3

g = g1 %u% g2 %u% g3
