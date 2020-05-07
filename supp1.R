library(reshape2)
library(ggplot2)
library(ggpubr)
library(BBmisc)
require(preprocessCore)
healthy = read.csv("healthy.csv")
disease = read.csv("disease.csv")
healthyforplot = healthy[, -1]
h = scale(healthyforplot) %>% as.data.frame()
hhh = disease[, -1]
h1 = scale(hhh) %>% as.data.frame()
h = cbind(healthyforplot, hhh)
h = as.matrix(h)
x = c(colnames(healthy)[-1], colnames(disease)[-1])
h = normalize.quantiles(h)
colnames(h) = x




g1 = ggplot(data = melt(h), aes(x=Var2, y= value)) + stat_boxplot(geom='errorbar')+ geom_boxplot(fill='#A4A4A4', color = "black") + theme(axis.text.x=element_blank(), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black")) + xlab("GSE65204 Samples")
g2 = ggplot(data = melt(h1), aes(x=variable, y= value)) + stat_boxplot(geom='errorbar')+ geom_boxplot(fill='#A4A4A4', color = "black")+ theme(axis.text.x=element_blank(), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black")) + xlab("Case")


plot = ggarrange(g1, g2, labels = c("A", "B"), ncol = 1, nrow = 2)
ggsave("fig1.tiff", plot = plot, dpi = 300)
