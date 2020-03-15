
##Usage of MetaIntegrator for Effect Size test
require(MetaIntegrator)
dataObj = cohortobj$originalData$GSE35571
dataObj = classFunction(dataObj, column = "disease status:ch1", diseaseTerms = "asthmatic")
dataObj$pheno$group = as.numeric(dataObj$class)

#Three datasets from MetaIntegrator
master = getGEOData(c( "GSE65204", "GSE64913", "GSE63142"))


#master = getGEOData(c( "GSE65204", "GSE64913", "GSE63142", "GSE40732", "GSE35571"))
data = master
##Cleaning up a weird dataset
gse65204 = read.csv("gse65204.csv")
classes = gse65204[8,]
classes = classes[, -1]
temp = unlist(classes)
for(i in 1:length(classes)){
	if(grepl("TRUE", temp[i]) == TRUE){
		classes[i] = 1
	}else{
		classes[i] = 0
	}
}
data$originalData$GSE65204$class = unlist(classes)


#Matching Classes with binaries
data$originalData$GSE64913 = classFunction(data$originalData$GSE64913, column = "diagnosis:ch1", diseaseTerms = c("Severe Asthmatic"))##OK....all these datasets are weird
#data$originalData$GSE35571 = classFunction(data$originalData$GSE35571, column = "disease status:ch1", diseaseTerms = c("asthmatic"))
data$originalData$GSE63142 = classFunction(data$originalData$GSE63142, column = "individual:ch1", diseaseTerms = c("severe asthmatic", "not severe asthmatic"))
dataObj = data$originalData$GSE63142
dataObj$pheno$group = as.numeric(dataObj$class)
data$originalData$GSE63142 = NULL
analysis = runMetaAnalysis(data, runLeaveOneOutAnalysis = F)
analysis = filterGenes(analysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.1)


##Plots for the meta analysis
violinPlot(analysis$filterResults$FDR0.05_es0_nStudies1_looaFALSE_hetero0, dataObj, labelColumn = 'group')
rocPlot(analysis$filterResults$FDR0.05_es0_nStudies1_looaFALSE_hetero0, dataObj)
pos = analysis$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$posGeneNames
neg = analysis$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$negGeneNames
deg = c(pos, neg)
results = analysis$metaAnalysis$pooledResults
sig = results[results$effectSizeFDR < 0.05,]
#Writing this file to a csv for later use
write.csv(sig, "degenes1.csv")


#Begin main script
#differential network creation
require(igraph)
require(reshape2)
require(ggplot2)
require(WGCNA)
require(flashClust)
require(org.Hs.eg.db)



#Loading in TF database
trr = read.csv("tfdb.csv")
tfs = unique(trr$targetID)
#Loading in DE genes from MetaIntegrator
m = read.csv("degenes1.csv", header = FALSE)
chars = as.character(m$V1)
#Mapping gene names to EENTREZ ID
omega = mapIds(org.Hs.eg.db, chars, 'ENTREZID', 'SYMBOL')
vec = as.numeric(omega)
#Removing missing values
vec = vec[!is.na(vec)]

#Loading in microarray data, taken from GEO
#A is for healthy, B for asthma
# GSE64913
datExprA1 = process(read.csv("healthy.csv"))
datExprB1 = process(read.csv("disease.csv"))
# GSE63142
datExprA2 = process(read.csv("healthy63.csv"))
datExprB2 = process(read.csv("disease63.csv"))
# GSE65204
datExprA3 = process(read.csv("healthy65.csv"))
datExprB3 = process(read.csv("disease65.csv"))


##WGCNA Gene Filtration Technique Applied

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


##Finding Thresholds for Gene Differential Co-expression networks

#Finding differential co-expression network
#GSE64913

adj1 = cor(t(datExprA1g))
adj2 = cor(t(datExprB1g))
adj = abs(adj2 - adj1)

#Thresholding the networks using cutoff identified in Figure 3
adj[adj < 1.3] = NA #lower if needed
df = melt(adj, na.rm = TRUE)
df = df[df$Var1 %in% vec,]
g1 = graph_from_data_frame(df)
denn1 = V(g1)$name

##GSE64913
adj1 = cor(t(datExprA2g))
adj2 = cor(t(datExprB2g))
adj = abs(adj2 - adj1)

#Thresholding the networks using cutoff identified in Figure 3
adj[adj < 0.9] = NA #lower if needed
df = melt(adj, na.rm = TRUE)
df = df[df$Var1 %in% vec,]
g2 = graph_from_data_frame(df)
denn2 = V(g2)$name


##GSE65142
adj1 = cor(t(datExprA3g))
adj2 = cor(t(datExprB3g))
adj = abs(adj2 - adj1)

#Thresholding the networks using cutoff identified in Figure 3
adj[adj < 0.9] = NA #lower if needed
df = melt(adj, na.rm = TRUE)
df = df[df$Var1 %in% vec,]
g3 = graph_from_data_frame(df)
denn3 = V(g3)$name

#Combining the three networks:

#individual gene lists (differentially expressed)
denn = denn1 %u% denn2 %u% denn3
#individual co-expression networks
g = g1 %u% g2 %u% g3


##Creation of Hive plots
k = dimnames(datExprA1g)[[1]]
#Finding only genes that are differentially expressed and/or a transcription factor
index = k %in% denn | k %in% tfs
datExprA1g = datExprA1g[index,]
datExprA2g = datExprA2g[index,] 
datExprA3g = datExprA3g[index,]
datExprB1g = datExprB1g[index,] 
datExprB2g = datExprB2g[index,] 
datExprB3g = datExprB3g[index,]


#Recreating differential network, setting insignificant values to 0.
cor1 = abs(cor(t(datExprA1g)) - cor(t(datExprB1g)))
cor2 = abs(cor(t(datExprA2g)) - cor(t(datExprB2g)))
cor3 = abs(cor(t(datExprA3g)) - cor(t(datExprB3g)))
adj = pmax(pmax(cor1, cor2), cor3)
adj[adj < 1.1] = 0
adj[adj > 0] = 1

#Finding interactions for GRN, combination of DE and tf interactions
trr = read.csv("trrust.csv")
dat = melt(adj)
trr = trr[trr$target %in% genes,]
temp = dat[dat$Var1 %in% tfs & dat$Var2 %in% trr$target, ]
temp = temp[temp$value == 1,]
temp1 = dat[dat$Var1 %in% genes & dat$Var2 %in% genes, ]
df = rbind(temp, temp1)
df = df[df$value == 1, ]

df$Var1 = df$Var1 %>% as.numeric
df$Var2 = df$Var2 %>% as.numeric

#Creation of hive plot object
hhp = edge2HPD(df, axis.cols = c("red", "orange", "green"))

deg = degree(g)
hhp$nodes$radius = 0
for(i in 1:length(hhp$nodes$radius)){
	label = hhp$nodes$lab[i]
	if(label %in% V(g)$name){
		hhp$nodes[i, ]$radius = deg[names(deg) == label] %>% as.numeric
	}
}

#Plotting hive plot.
hhp$nodes$radius = hhp$nodes$radius %>% jitter()
hhp = manipAxis(hhp, 'norm')
hhp$nodes$color = "white"
hhp$nodes$size = 0.5
hhp$edges$weight = 0.5
hhp = manipAxis(hhp, 'norm')
