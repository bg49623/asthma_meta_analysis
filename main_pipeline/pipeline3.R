### Pipeline 1.3.
# Things done here:
# Module Preservation across healthy and asthma clusters using WGCNA framework
# Fisher's test to determine asthma-unique modules
# Plotting of heatmap to identify asthma-unique modules
# GO term analysis using goana (limma) function


##Load in Asthmatic Mains


B1 = readRDS("treeObjB1.rds")
B2 = readRDS("treeObjB2.rds")
B3 = readRDS("treeObjB3.rds")


#Asthma Base Network here.
modulesA1 = B1[[1]]
mappings = data.frame(B1[[1]], B1[[2]])
colnames(mappings) = c("color", "gene")
multiExpr = list(A1=list(data=t(datExprB1g)),A2=list(data=t(datExprA1g)))
multiColor = list(A1 = modulesA1)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=50,maxGoldModuleSize=100,maxModuleSize=400)
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
ranks = mp$preservation$observed$ref.A1$inColumnsAlsoPresentIn.A2
ranking = stats[order(-stats[,2]),c(1:2)]
lowpres = ranking[ranking$Zsummary.pres < 10, ]
preservedA = row.names(lowpres)
preservedA = preservedA[!preservedA %in% c("gold", "grey")]
prescolorsA = mappings[mappings$color %in% preservedA,]
prescolorsA$color = as.factor(prescolorsA$color)
goodmodsA = split(prescolorsA, f= prescolorsA$color)
saveRDS(goodmodsA, "asthmamods1.rds")


modulesA1 = B2[[1]]
mappings = data.frame(B2[[1]], B2[[2]])
colnames(mappings) = c("color", "gene")
multiExpr = list(A1=list(data=t(datExprB2g)),A2=list(data=t(datExprA2g)))
multiColor = list(A1 = modulesA1)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=50,maxGoldModuleSize=100,maxModuleSize=400)
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
ranks = mp$preservation$observed$ref.A1$inColumnsAlsoPresentIn.A2
ranking = stats[order(-stats[,2]),c(1:2)]
lowpres = ranking[ranking$Zsummary.pres < 10, ]
preservedA = row.names(lowpres)
preservedA = preservedA[!preservedA %in% c("gold", "grey")]
prescolorsA = mappings[mappings$color %in% preservedA,]
prescolorsA$color = as.factor(prescolorsA$color)
goodmodsA = split(prescolorsA, f= prescolorsA$color)
saveRDS(goodmodsA, "asthmamods2.rds")


modulesA1 = B3[[1]]
mappings = data.frame(B3[[1]], B3[[2]])
colnames(mappings) = c("color", "gene")
multiExpr = list(A1=list(data=t(datExprB3g)),A2=list(data=t(datExprA3g)))
multiColor = list(A1 = modulesA1)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=50,maxGoldModuleSize=100,maxModuleSize=400)
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
ranks = mp$preservation$observed$ref.A1$inColumnsAlsoPresentIn.A2
ranking = stats[order(-stats[,2]),c(1:2)]
lowpres = ranking[ranking$Zsummary.pres < 10, ]
preservedA = row.names(lowpres)
preservedA = preservedA[!preservedA %in% c("gold", "grey")]
prescolorsA = mappings[mappings$color %in% preservedA,]
prescolorsA$color = as.factor(prescolorsA$color)
goodmodsA = split(prescolorsA, f= prescolorsA$color)
saveRDS(goodmodsA, "asthmamods3.rds")

library(stringr)
a = readRDS("asthmamods1.rds")
b = readRDS("asthmamods2.rds")
c = readRDS("asthmamods3.rds")

a = Filter(function(x) {length(x$gene) > 0}, a)
names(a) = str_replace_all(paste(names(a), "A"), " ", "")
for(i in 1:length(names(a))){
  str = names(a)[i]
  a[[i]]$color = str
}
b = Filter(function(x) {length(x$gene) > 0}, b)
names(b) = str_replace_all(paste(names(b), "B"), " ", "")
for(i in 1:length(names(b))){
  str = names(b)[i]
  b[[i]]$color = str
}
c = Filter(function(x) {length(x$gene) > 0}, c)
names(c) = str_replace_all(paste(names(c), "C"), " ", "")
for(i in 1:length(names(c))){
  str = names(c)[i]
  c[[i]]$color = str
}
# filtering for the nonempty elements of the list. 

main = c(a, b, c)
answers = c()
for(i in main){
  for(j in main){
    
    pval = fisher1(i,j)
    answers = c(answers,pval)
  }
}

lena = length(main)
mat1 = matrix(answers, nrow = lena)
rownames(mat1) = names(main)
colnames(mat1) = names(main)
df = melt(mat1)
df$val = str_sub(df$Var1, start = -1)
df$val = factor(df$val)
g1 = ggplot(data = df,  aes(Var2, Var1, fill = value))+ geom_tile(color = "black") + scale_fill_gradient2(low="black", high="blue", guide="colorbar") + coord_equal() + theme_classic()+
  theme(axis.ticks = element_blank(), axis.text.x=element_text(angle = 90, hjust = 1),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.title.x = element_text(color = "black", size = 6, face = "bold"),
        axis.text.y = element_text(color = "black", size = 6, face = "bold")) + facet_wrap(~val) + xlab("Modules (x3)")


total_genes = 12806


setsA = list()
answers = c()
targets1 = c()
targets2 = c()
counter = 1
for(i in a){
  for(j in b){
    for(k in c){
      pval = fisher(i, j)
      pval1 = fisher(i, k)
      if(pval < 0.05 & pval1 < 0.05){
        l = unique(c(i$gene, j$gene, k$gene))
        setsA[[counter]] = l
        counter = counter + 1
        answers = c(answers, as.character(i$color[1]))
        targets1 = c(targets1, as.character(j$color[1]))
        targets2 = c(targets2, as.character(k$color[1]))
      }
    }
  }
}

df1 = data.frame(answers, targets1)
df2 = data.frame(answers, targets2)
colnames(df2) = c("answers", "targets1")
dfA = rbind(df1, df2)

for(i in 1:length(setsA)){
  temp = commonProbesA[setsA[[i]]]
  setsA[[i]] = temp
}

setsB = list()
answers = c()
targets1 = c()
targets2 = c()
counter = 1
for(i in b){
  for(j in a){
    for(k in c){
      pval = fisher(i, j)
      pval1 = fisher(i, k)
      if(pval < 0.05 & pval1 < 0.05){
        l = unique(c(i$gene, j$gene, k$gene))
        setsB[[counter]] = l
        counter = counter + 1
        answers = c(answers, as.character(i$color[1]))
        targets1 = c(targets1, as.character(j$color[1]))
        targets2 = c(targets2, as.character(k$color[1]))
      }
    }
  }
}

for(i in 1:length(setsB)){
  temp = commonProbesA[setsB[[i]]]
  setsB[[i]] = temp
}

df1 = data.frame(answers, targets1)
df2 = data.frame(answers, targets2)
colnames(df2) = c("answers", "targets1")
dfB = rbind(df1, df2)



setsC = list()
answers = c()
targets1 = c()
targets2 = c()
counter = 1
for(i in c){
  for(j in b){
    for(k in a){
      pval = fisher(i, j)
      pval1 = fisher(i, k)
      if(pval < 0.05 & pval1 < 0.05){
        l = unique(c(i$gene, j$gene, k$gene))
        setsC[[counter]] = l
        counter = counter + 1
        answers = c(answers, as.character(i$color[1]))
        targets1 = c(targets1, as.character(j$color[1]))
        targets2 = c(targets2, as.character(k$color[1]))
      }
    }
  }
}
for(i in 1:length(setsC)){
  temp = commonProbesA[setsC[[i]]]
  setsC[[i]] = temp
}

df1 = data.frame(answers, targets1)
df2 = data.frame(answers, targets2)
colnames(df2) = c("answers", "targets1")
dfB = rbind(df1, df2)

df1 = data.frame(answers, targets1)
df2 = data.frame(answers, targets2)
colnames(df2) = c("answers", "targets1")
dfC = rbind(df1, df2) %>% unique()


df = rbind(dfA, dfB, dfC)

fisher = function(i, j){
  intersection = length(intersect(i$gene,j$gene))
  if(intersection == 0){
    return(1)
  }
  topright = length(i$gene) - intersection
  botleft = length(j$gene) - intersection
  botright = total_genes - topright - botleft - intersection
  mat = matrix(c(intersection, topright, botleft, botright), ncol = 2)
  pval = fisher.test(mat)$p.value
  return(pval)
}
fisher1 = function(i, j){
  intersection = length(intersect(i$gene,j$gene))
  if(intersection == 0){
    return(1)
  }
  topright = length(i$gene) - intersection
  botleft = length(j$gene) - intersection
  botright = total_genes - topright - botleft - intersection
  mat = matrix(c(intersection, topright, botleft, botright), ncol = 2)
  pval = fisher.test(mat)$p.value
  if(pval < 0.05){
    return(0)
  }
  return(1)
}


##Reduce the answers.
answers1 = unique(answers1)
answers2 = unique(answers2)
answers3 = unique(answers3)

for(i in a){
  if(length(i$gene)> 0){
    for(j in b){
      if(length(j$gene)> 0){
        # Need to get counts for Fisher's exact test
        # We construct a 2x2 matrix, where the top left value is the intersction of the two modules
        # The top right and bottom left are the unique values to each module. 
        # The bottom right is everything that is left
        intersection = length(intersect(i$gene,j$gene))
        topright = length(i$gene) - intersection
        botleft = length(j$gene) - intersection
        botright = total_genes - topright - botleft - intersection
        mat = matrix(c(intersection, topright, botleft, botright), ncol = 2)
        pval = fisher.test(mat)$p.value
      }
    }
    for(k in c){
      if(length(k$gene)> 0){
        # Need to get counts for Fisher's exact test
        # We construct a 2x2 matrix, where the top left value is the intersction of the two modules
        # The top right and bottom left are the unique values to each module. 
        # The bottom right is everything that is left
        intersection = length(intersect(i$gene,j$gene))
        topright = length(i$gene) - intersection
        botleft = length(j$gene) - intersection
        botright = total_genes - topright - botleft - intersection
        mat = matrix(c(intersection, topright, botleft, botright), ncol = 2)
        pval1 = fisher.test(mat)$p.value
      }
    }

  }
}


for(i in a){
  print(i)
}
pvals = matrix(pval, ncol = 5)
print(pvals)



## Base: Singhania et al.
answers = c()
answers1 = c()

for(i in a){
  for(j in b){
    pval = fisher(i, j)
    answers = c(answers, pval)
    
  }
}
for(i in a){
  for(k in c){
    pval = fisher(i, k)
    answers1 = c(answers1, pval)
    
  }
}
lena = length(a)
mat1 = matrix(answers, nrow = lena)
rownames(mat1) = names(a)
mat2 = matrix(answers1, nrow = lena)
rownames(mat2) = names(a)
g1 = ggplot(data = melt(mat1),  aes(Var2, Var1, fill = value))+ geom_tile(color = "black") + scale_fill_gradient2(low="black", high="blue", guide="colorbar") + coord_equal() + theme_classic()+
  theme(axis.ticks = element_blank(), axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, face = "bold")) + xlab("Wenzel et al.")
g2 = ggplot(data = melt(mat2),  aes(Var2, Var1, fill = value))+ geom_tile(color = "black") + scale_fill_gradient2(low="black", high="blue", guide="colorbar") + coord_equal() + theme_classic()+
  theme(axis.ticks = element_blank(), axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, face = "bold")) + xlab("Yang et al.")


## Base: Wenzel et al.
answers = c()
answers1 = c()

for(i in b){
  for(j in a){
    pval = fisher(i, j)
    answers = c(answers, pval)
    
  }
}
for(i in b){
  for(k in c){
    pval = fisher(i, k)
    answers1 = c(answers1, pval)
    
  }
}
lena = length(b)
mat1 = matrix(answers, nrow = lena)
rownames(mat1) = names(b)
mat2 = matrix(answers1, nrow = lena)
rownames(mat2) = names(b)
g1 = ggplot(data = melt(mat1),  aes(Var2, Var1, fill = value))+ geom_tile(color = "black") + scale_fill_gradient2(low="black", high="blue", guide="colorbar") + coord_equal() + theme_classic()+
  theme(axis.ticks = element_blank(), axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, face = "bold")) + xlab("Singhania et al.")
g2 = ggplot(data = melt(mat2),  aes(Var2, Var1, fill = value))+ geom_tile(color = "black") + scale_fill_gradient2(low="black", high="blue", guide="colorbar") + coord_equal() + theme_classic()+
  theme(axis.ticks = element_blank(), axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, face = "bold")) + xlab("Yang et al.")

## Base: Yang et al.
answers = c()
answers1 = c()

for(i in c){
  for(j in a){
    pval = fisher(i, j)
    answers = c(answers, pval)
    
  }
}
for(i in c){
  for(k in b){
    pval = fisher(i, k)
    answers1 = c(answers1, pval)
    
  }
}
lena = length(c)
mat1 = matrix(answers, nrow = lena)
rownames(mat1) = names(c)
mat2 = matrix(answers1, nrow = lena)
rownames(mat2) = names(c)
g1 = ggplot(data = melt(mat1),  aes(Var2, Var1, fill = value))+ geom_tile(color = "black") + scale_fill_gradient2(low="black", high="blue", guide="colorbar") + coord_equal() + theme_classic()+
  theme(axis.ticks = element_blank(), axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, face = "bold")) + xlab("Singhania et al.")
g2 = ggplot(data = melt(mat2),  aes(Var2, Var1, fill = value))+ geom_tile(color = "black") + scale_fill_gradient2(low="black", high="blue", guide="colorbar") + coord_equal() + theme_classic()+
  theme(axis.ticks = element_blank(), axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
        axis.text.y = element_text(color = "black", size = 8, face = "bold")) + xlab("Yang et al.")
