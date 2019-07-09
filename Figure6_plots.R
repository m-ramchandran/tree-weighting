#First, run singledataset_examples.R

#(a) GSE25066: Binary outcome (dmfs status)
plot_means.z <- data.frame(Method = c("Unweighted", "Weighting Forests", "Weighting Trees", "Combined"), 
                           Means = cbind(rep_test$means_tree[c(11, 13)], 
                                         rep_test$means_forest[c(2, 12)])[c(3, 4,1,2)],
                           Sds = cbind(rep_test$sds_tree[c(11, 13)], 
                                       rep_test$sds_forest[c(2, 12)])[c(3,4,1,2)])

p2 <- ggplot(plot_means.z, aes(x=factor(Method, levels = c("Unweighted", "Weighting Forests", "Weighting Trees", "Combined")), y = Means, fill = Method)) + geom_bar(stat = "identity", fill = c("#F4A261","#254653", "#2A9D8F", "purple")) +
  theme_classic() + geom_errorbar(aes(ymin=Means-Sds/10, ymax=Means+Sds/10), width=.2) + xlab("Method") + 
  ylab("Percent change in average Log Loss from Merged") + labs(title = "GSE25066: Binary outcome (dmfs status)") + scale_y_continuous(breaks = seq(-100, 0, by=10), limits=c(-100,0)) +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)))

show(p2)


#(b) GSE25066: Continuous outcome (days to dmfs)
plot_means.z1 <- data.frame(Method = c("Unweighted", "Weighting Forests", "Weighting Trees", "Combined"), 
                            Means = cbind(rep_test1$means_tree[c(11, 13)], 
                                          rep_test1$means_forest[c(2, 12)])[c(3, 4,1,2)],
                            Sds = cbind(rep_test1$sds_tree[c(11, 13)], 
                                        rep_test1$sds_forest[c(2, 12)])[c(3,4,1,2)])

p3 <- ggplot(plot_means.z1, aes(x=factor(Method, levels = c("Unweighted", "Weighting Forests", "Weighting Trees", "Combined")), y = Means, fill = Method)) + geom_bar(stat = "identity", fill = c("purple","#F4A261", "#254653", "#2A9D8F")) +
  theme_classic() + geom_errorbar(aes(ymin=Means-Sds/10, ymax=Means+Sds/10), width=.2) + xlab("Method") +
  ylab("Percent change in average RMSE from Merged") + labs(title = "GSE25066: Continuous outcome (days to dmfs)")  +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)))

show(p3)


#(c) Metabric: Continuous outcome (days to death)
plot_means.z2 <- data.frame(Method = c("Unweighted", "Weighting Forests", "Weighting Trees", "Combined"), 
                            Means = cbind(rep_test4$means_tree[c(11, 13)], 
                                          rep_test4$means_forest[c(2, 12)])[c(3,4,1,2)],
                            Sds = cbind(rep_test4$sds_tree[c(11, 13)], 
                                        rep_test4$sds_forest[c(2, 12)])[c(3,4,1,2)])

library(wesanderson)

p4 <- ggplot(plot_means.z2, aes(x=factor(Method, levels = c("Unweighted", "Weighting Forests", "Weighting Trees", "Combined")), y = Means, fill = Method)) + geom_bar(stat = "identity") + scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) +#, fill = c("#F4A261","purple", "#2A9D8F", "#254653")) +
  theme_classic() + geom_errorbar(aes(ymin=Means-Sds/10, ymax=Means+Sds/10), width=.2) + theme(legend.position="none") + xlab("Method") +
  ylab("Percent change in average RMSE from Merged") + labs(title = "Metabric: Continuous outcome (days to death)")  +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)))

show(p4)


#Create Figure 6
library(gridExtra)
library(ggpubr)
ggarrange(p2, p3, p4, labels = c("A", "B", "C"), ncol=3, nrow=1)




