require(MetaIntegrator)
dataObj = cohortobj$originalData$GSE35571
dataObj = classFunction(dataObj, column = "disease status:ch1", diseaseTerms = "asthmatic")
dataObj$pheno$group = as.numeric(dataObj$class)

#Three datasets from MetaIntegrator
master = getGEOData(c( "GSE65204", "GSE64913", "GSE63142"))


#master = getGEOData(c( "GSE65204", "GSE64913", "GSE63142", "GSE40732", "GSE35571"))
data = master
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