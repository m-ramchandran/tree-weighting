#First, run main_simulation.R

#vector of interaction strength coefficients
inter_ls <- seq(2, 12, .5)

#vector of heterogeneity 'bad' level coefficients
het_ls <- seq(0, 10, .5)

#Figure 2
#(a)
#reading in files
means_c1s1 <- read.csv("means_c1_tree.csv")
sds_c1s1 <- read.csv("sds_c1_tree.csv")
means_c2s1 <- read.csv("means_c1_forest.csv")
sds_c2s1 <- read.csv("sds_c1_forest.csv")

#Creating the plotting data frames
cont_s1 <- as.matrix(cbind(means_c1s1[,c(2, 13)], means_c2s1[,c(1, 13)])[,c(3,1,2,4)])
data_c1 <- data.frame(cbind(inter_ls, cont_s1))
df.melted1 <- melt(data_c1, id = "inter_ls")
se_c1.Tree <- data.frame(cbind(inter_ls, means_c1s1$Stack_ridge - 1.96*sds_c1s1$Stack_ridge/(sqrt(100)),
                               means_c1s1$Stack_ridge + 1.96*sds_c1s1$Stack_ridge/(sqrt(100))))
se_c1.Forest <- data.frame(cbind(inter_ls, means_c2s1$Stack_ridge - 1.96*sds_c2s1$Stack_ridge/(sqrt(100)),
                                 means_c2s1$Stack_ridge + 1.96*sds_c2s1$Stack_ridge/(sqrt(100))))
se_c1.Unweighted <- data.frame(cbind(inter_ls, means_c1s1$Unweighted - 1.96*sds_c1s1$Unweighted/(sqrt(100)),
                                     means_c1s1$Unweighted + 1.96*sds_c1s1$Unweighted/(sqrt(100))))
se_c1.Merged <- data.frame(cbind(inter_ls, means_c2s1$Merged - 1.96*sds_c2s1$Merged/(sqrt(100)),
                                 means_c2s1$Merged + 1.96*sds_c2s1$Merged/(sqrt(100))))
#Plot
q1 <- ggplot(data = df.melted1, aes(x = inter_ls, y = value, color = variable)) +
  geom_point() + stat_smooth(se = F) + 
  geom_ribbon(data=se_c1.Tree,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) + 
  geom_ribbon(data=se_c1.Forest,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#254653",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_c1.Unweighted,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#F4A261",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_c1.Merged,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#E76E51",alpha=0.2, inherit.aes = FALSE) +
  theme_classic() + xlab("Interaction strength") + ylab("Average RMSE") + 
  labs(title = "2 interactions in training, 2 in testing") +
  theme(legend.position="none") + ylim(4, 8) + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5))) +
  scale_color_manual(values=c("#E76E51","#F4A261","#2A9D8F","#254653"), labels = c("Merged", "Unweighted", "Weighting Trees", "Weighting Forests"), name="Method")
show(q1)

#(b)
#reading in files
means_c1s2 <- read.csv("~/Desktop/case2/means_c2_tree.csv")
sds_c1s2 <- read.csv("~/Desktop/case2/sds_c2_tree.csv")
means_c2s2 <- read.csv("~/Desktop/case2/means_c2_forest.csv")
sds_c2s2 <- read.csv("~/Desktop/case2/sds_c2_forest.csv")

#Creating the plotting data frames
cont_s2 <- as.matrix(cbind(means_c1s2[,c(2, 13)], means_c2s2[,c(1, 13)])[,c(3,1,2,4)])
data_c2 <- data.frame(cbind(inter_ls, cont_s2))
df.melted2 <- melt(data_c2, id = "inter_ls")
se_c2.Tree <- data.frame(cbind(inter_ls, means_c1s2$Stack_ridge - 1.96*sds_c1s2$Stack_ridge/(sqrt(100)),
                               means_c1s2$Stack_ridge + 1.96*sds_c1s2$Stack_ridge/(sqrt(100))))
se_c2.Forest <- data.frame(cbind(inter_ls, means_c2s2$Stack_ridge - 1.96*sds_c2s2$Stack_ridge/(sqrt(100)),
                                 means_c2s2$Stack_ridge + 1.96*sds_c2s2$Stack_ridge/(sqrt(100))))
se_c2.Unweighted <- data.frame(cbind(inter_ls, means_c1s2$Unweighted - 1.96*sds_c1s2$Unweighted/(sqrt(100)),
                                     means_c1s2$Unweighted + 1.96*sds_c1s2$Unweighted/(sqrt(100))))
se_c2.Merged <- data.frame(cbind(inter_ls, means_c2s2$Merged - 1.96*sds_c2s2$Merged/(sqrt(100)),
                                 means_c2s2$Merged + 1.96*sds_c2s2$Merged/(sqrt(100))))
#Plot
q2<- ggplot(data = df.melted2, aes(x = inter_ls, y = value, color = variable)) +
  geom_point() + stat_smooth(se = F) + 
  geom_ribbon(data=se_c2.Tree,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) + 
  geom_ribbon(data=se_c2.Forest,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#254653",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_c2.Unweighted,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#F4A261",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_c2.Merged,aes(x=inter_ls,ymin=V2,ymax=V3),fill="#E76E51",alpha=0.2, inherit.aes = FALSE) +
  theme_classic() + xlab("Interaction strength") + ylab("Average RMSE") + 
  labs(title = "6 interactions in training, 2 in testing") + 
  theme(legend.position="none") + ylim(4, 8) + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5))) +
  scale_color_manual(values=c("#E76E51","#F4A261","#2A9D8F","#254653"), labels = c("Merged", "Unweighted", "Weighting Trees", "Weighting Forests"), name="Method")
show(q2)

#(c) 
#reading in files
means_c1s3 <- read.csv("~/Desktop/case3/means_c3_tree.csv")
sds_c1s3 <- read.csv("~/Desktop/case3/sds_c3_tree.csv")
means_c2s3 <- read.csv("~/Desktop/case3/means_c3_forest.csv")
sds_c2s3 <- read.csv("~/Desktop/case3/sds_c3_forest.csv")

#Creating the plotting data frames
cont_s3 <- as.matrix(cbind(means_c1s3[,c(2, 13)], means_c2s3[,c(1, 13)])[,c(3,1,2,4)])
data_c3 <- data.frame(cbind(het_ls, cont_s3))
df.melted3 <- melt(data_c3, id = "het_ls")
se_c3.Tree <- data.frame(cbind(het_ls, means_c1s3$Stack_ridge - 1.96*sds_c1s3$Stack_ridge/(sqrt(100)),
                               means_c1s3$Stack_ridge + 1.96*sds_c1s3$Stack_ridge/(sqrt(100))))
se_c3.Forest <- data.frame(cbind(het_ls, means_c2s3$Stack_ridge - 1.96*sds_c2s3$Stack_ridge/(sqrt(100)),
                                 means_c2s3$Stack_ridge + 1.96*sds_c2s3$Stack_ridge/(sqrt(100))))
se_c3.Unweighted <- data.frame(cbind(het_ls, means_c1s3$Unweighted - 1.96*sds_c1s3$Unweighted/(sqrt(100)),
                                     means_c1s3$Unweighted + 1.96*sds_c1s3$Unweighted/(sqrt(100))))
se_c3.Merged <- data.frame(cbind(het_ls, means_c2s3$Merged - 1.96*sds_c2s3$Merged/(sqrt(100)),
                                 means_c2s3$Merged + 1.96*sds_c2s3$Merged/(sqrt(100))))
#Plot
q3 <- ggplot(data = df.melted3, aes(x = het_ls, y = value, color = variable)) +
  geom_point() + stat_smooth(se = F) + 
  geom_ribbon(data=se_c3.Tree,aes(x=het_ls,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) + 
  geom_ribbon(data=se_c3.Forest,aes(x=het_ls,ymin=V2,ymax=V3),fill="#254653",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_c3.Unweighted,aes(x=het_ls,ymin=V2,ymax=V3),fill="#F4A261",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_c3.Merged,aes(x=het_ls,ymin=V2,ymax=V3),fill="#E76E51",alpha=0.2, inherit.aes = FALSE) +
  theme_classic() + xlab("Heterogeneity Level") + ylab("Average RMSE") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5))) +
  labs(title = "No datasets with interactions") + 
  scale_color_manual(values=c("#E76E51","#F4A261","#2A9D8F","#254653"), labels = c("Merged", "Unweighted", "Weighting Trees", "Weighting Forests"), name="Method")
show(q3)

#Create Figure 2:
library(gridExtra)
library(ggpubr)
ggarrange(q1,q2, q3, labels = c("A", "B", "C"), ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

