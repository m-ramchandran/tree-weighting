#First, run tcga_combined.R

#vector of heterogeneity 'bad' level coefficients
het_ls <- seq(0, 10, .5)

#reading in files
means_c1s4 <- read.csv("~/Desktop/tcga_moretrees/means_c1s4_tree.csv")
sds_c1s4 <- read.csv("~/Desktop/tcga_moretrees/sds_c1s4_tree.csv")
means_c2s4 <- read.csv("~/Desktop/tcga_moretrees/means_c1s4_forest.csv")
sds_c2s4 <- read.csv("~/Desktop/tcga_moretrees/sds_c1s4_forest.csv")


#Creating the plotting data frames
cont_tcga <- as.matrix(cbind(means_c1s4$Unweighted,means_c1s4$Stack_ridge, means_c2s4$Stack_ridge, means_c1s4$Stack_ridge_combined))
data_tcga <- data.frame(cbind(het_ls, cont_tcga))
df.melted4 <- melt(data_tcga, id = "het_ls")
se_tcga.Tree <- data.frame(cbind(het_ls, means_c1s4$Stack_ridge - 1.96*sds_c1s4$Stack_ridge/(sqrt(100)),
                                 means_c1s4$Stack_ridge + 1.96*sds_c1s4$Stack_ridge/(sqrt(100))))
se_tcga.Forest <- data.frame(cbind(het_ls, means_c2s4$Stack_ridge - 1.96*sds_c2s4$Stack_ridge/(sqrt(100)),
                                   means_c2s4$Stack_ridge + 1.96*sds_c2s4$Stack_ridge/(sqrt(100))))
se_tcga.Unweighted <- data.frame(cbind(het_ls, means_c1s4$Unweighted - 1.96*sds_c1s4$Unweighted/(sqrt(100)),
                                       means_c1s4$Unweighted + 1.96*sds_c1s4$Unweighted/(sqrt(100))))
se_tcga.Combined <- data.frame(cbind(het_ls, means_c1s4$Stack_ridge_combined - 1.96*sds_c1s4$Stack_ridge_combined/(sqrt(100)),
                                     means_c1s4$Stack_ridge_combined + 1.96*sds_c1s4$Stack_ridge_combined/(sqrt(100))))

#Plot
q6 <- ggplot(data = df.melted4, aes(x = het_ls, y = value, color = variable)) +
  geom_point() + stat_smooth(se = F) + 
  geom_ribbon(data=se_tcga.Tree,aes(x=het_ls,ymin=V2,ymax=V3),fill="darkcyan",alpha=0.2, inherit.aes = FALSE) + 
  geom_ribbon(data=se_tcga.Forest,aes(x=het_ls,ymin=V2,ymax=V3),fill="midnightblue",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_tcga.Unweighted,aes(x=het_ls,ymin=V2,ymax=V3),fill="sienna2",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_tcga.Combined,aes(x=het_ls,ymin=V2,ymax=V3),fill="darkorchid1",alpha=0.2, inherit.aes = FALSE) +
  theme_classic() + xlab("Heterogeneity Level") + ylab("Percent change in average RMSE from Merged") +
  theme(legend.position="bottom") +
  labs(title = "TCGA study") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5))) +
  scale_color_manual(values=c("sienna2","darkcyan", "midnightblue", "darkorchid1"), labels = c("Unweighted", "Weighting Trees", "Weighting Forests", "Combined"), name="Method")
#scale_color_manual(values=c("red","darkgreen","turquoise","purple"), labels = c("Merged", "Unweighted", "Weighting Trees", "Weighting Forests"), name="Method")
show(q6)
