### Pipeline 1.0
# Using MetaIntegrator on our three datasets to determine list of Differentially expressed genes. 
# Testing these genes with binary classification tests (ROC) for validation


require(MetaIntegrator)
#master = getGEOData(c( "GSE65204", "GSE64913", "GSE35571", "GSE63142"))
dataObj = cohortobj$originalData$GSE35571
dataObj = classFunction(dataObj, column = "disease status:ch1", diseaseTerms = "asthmatic")
dataObj$pheno$group = as.numeric(dataObj$class)
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





#data$originalData$GSE40732 = classFunction(data$originalData$GSE40732, column = "asthma:ch1", diseaseTerms = c("TRUE"))
data$originalData$GSE64913 = classFunction(data$originalData$GSE64913, column = "diagnosis:ch1", diseaseTerms = c("Severe Asthmatic"))##OK....all these datasets are weird
#data$originalData$GSE35571 = classFunction(data$originalData$GSE35571, column = "disease status:ch1", diseaseTerms = c("asthmatic"))
data$originalData$GSE63142 = classFunction(data$originalData$GSE63142, column = "individual:ch1", diseaseTerms = c("severe asthmatic", "not severe asthmatic"))
dataObj = data$originalData$GSE63142
dataObj$pheno$group = as.numeric(dataObj$class)
data$originalData$GSE63142 = NULL
analysis = runMetaAnalysis(data, runLeaveOneOutAnalysis = F)
analysis = filterGenes(analysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.1)

violinPlot(analysis$filterResults$FDR0.05_es0_nStudies1_looaFALSE_hetero0, dataObj, labelColumn = 'group')
rocPlot(analysis$filterResults$FDR0.05_es0_nStudies1_looaFALSE_hetero0, dataObj)
pos = analysis$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$posGeneNames
neg = analysis$filterResults$FDR0.1_es0_nStudies1_looaFALSE_hetero0$negGeneNames
deg = c(pos, neg)
results = analysis$metaAnalysis$pooledResults
sig = results[results$effectSizeFDR < 0.05,]




write.csv(sig, "degenes1.csv")

m = read.csv("degenes1.csv", header = T)
library(org.Hs.eg.db)
chars = as.character(m$x)
omega = mapIds(org.Hs.eg.db, chars, 'ENTREZID', 'SYMBOL')
vec = as.numeric(omega)
vec = vec[!is.na(vec)]
trr = read.csv("tfdb.csv")
##Need to Map gene codes to ENTREZ
trrust = read.csv("trrust.csv")
l1 = as.character(trrust$A)
l2 = as.character(trrust$B)
omega1 = as.numeric(mapIds(org.Hs.eg.db, l1, 'ENTREZID', 'SYMBOL'))
omega2 = as.numeric(mapIds(org.Hs.eg.db, l2, 'ENTREZID', 'SYMBOL'))
trrust$A = omega1
trrust$B = omega2
trrust = na.omit(trrust)
trr = na.omit(trr)
##Where A is the TF and B is the target gene.
##These two are constant for both of the datasets: only need to run their respective scripts once. 

##Setting up correlation matrix
library(reshape2)
library(dplyr)
##Normalizing Gene expression Matrix
healthy = read.csv("healthy.csv")
vals = (healthy[, -1])
vals = scale(vals)
mapping = as.numeric(healthy[, 1])
m = data.matrix(vals)
mt = t(m)
cors = cor(mt)
index = as.numeric()
mappings = as.numeric()
##We only care about mapping to certain genes. This constricts it to only the DE genes
for(i in 1:length(as.numeric(mapping))){
	if(mapping[i] %in% vec){
		index = c(index, i)
		mappings = c(mappings, mapping[i])
	}
}
##Gives edgelist with all correlations to DE genes
cors = cors[index, ]
colnames(cors) = mapping
rownames(cors) = mappings
#cors[cors < 0] = NA
#Total Edgelist
edge = melt(cors, na.rm = TRUE)

#External TF database for checking
##In edge, we have Var2 -> Var1, where length(unique(edge$Var1)) has the 685 genes of interest
##Performing this in intermediate steps for memory reasons. 
##First truncate TRRUST list to only include our DE genes as targets
require(dplyr)
colnames(trr) = c("A", "B")
tfactors = unique(trr$B)
edge1 = edge[(edge$Var1 %in% tfactors ),]
temp = unique(edge1$Var2)
edge2 = edge[(edge$Var1 %in% temp & edge$Var2 %in% temp),]
edge2 = edge2[edge2$value > 0.9,]
hsmall = edge1
healthy_edgelist = rbind(edge1, edge2)
colnames(healthy_edgelist) = c("V1", "V2", "weight")
write.csv(x, "edgelist1.csv", row.names = F)
#Where A in x are the transcription factors, and B represents the targets of the TFs

##Now the other dataset
dis = read.csv("disease.csv")
vals = (dis[, -1])
vals = scale(vals)
mapping = as.numeric(dis[, 1])
m = data.matrix(vals)
mt = t(m)
cors = cor(mt)
index = as.numeric()
mappings = as.numeric()
##We only care about mapping to certain genes. This constricts it to only the DE genes
for(i in 1:length(as.numeric(mapping))){
	if(mapping[i] %in% vec){
		index = c(index, i)
		mappings = c(mappings, mapping[i])
	}
}

cors = cors[index, ]
colnames(cors) = mapping
rownames(cors) = mappings
#cors[cors < 0] = NA
edge = melt(cors, na.rm = TRUE)
colnames(trr) = c("A", "B")
tfactors = unique(trr$B)
edge1 = edge[(edge$Var1 %in% tfactors ),]
temp = unique(edge1$Var2)
edge2 = edge[(edge$Var1 %in% temp & edge$Var2 %in% temp),]
edge2 = edge2[edge2$value > 0.9,]
dsmall = edge1

disease_edgelist = rbind(edge1, edge2)
colnames(disease_edgelist) = c("V1", "V2", "weight")
write.csv(x, "edgelist2.csv", row.names = F)

library(igraph)
x = healthy_edgelist
w = x[(x$weight > 0.85),]
total  = unique(c(w$V1, w$V2))
hnet = graph_from_data_frame(d= w[, 1:2], vertices = total, directed = T)
hnet = simplify(hnet, remove.multiple = T, remove.loops = T)

##Identify the TF genes, DE genes, and Other genes in each network. 
##First, find TF's, these are the sources in our TFDB.

trr = read.csv("tfdb.csv")
tfs = trr$targetID
##List of TFs: this is constant for both the networks, so we only have to do it once. 
vertices = vertex_attr(hnet)$name %>% as.numeric()
htflist = vertices[vertices %in% tfs] %>% as.numeric()
htarglist = vertices[vertices %in% vec & !(vertices %in% htflist)]

hotherlist = vertices[!(vertices %in% c(htflist,htarglist))]
y = disease_edgelist
wd = y[(y$weight > 0.75), ]
total  = unique(c(wd$V1, wd$V2))
dnet = graph_from_data_frame(d= wd[, 1:2], vertices = total, directed = T)
dnet = simplify(dnet, remove.multiple = T, remove.loops = T)
vertices1 = vertex_attr(dnet)$name %>% as.numeric()
dtflist = vertices1[vertices1 %in% tfs] %>% as.numeric()
dtarglist = vertices1[vertices1 %in% vec & !(vertices1 %in% dtflist)]
dotherlist = vertices1[!(vertices1 %in% c(dtflist,dtarglist))]
