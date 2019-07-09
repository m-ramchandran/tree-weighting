#(a) % variation in the outcome explained by interaction terms
#First, run variance_inter_sim.R

#plot
inter_ls <- seq(2, 12, .5)

#data 
vi_data <- as.data.frame(cbind(inter_ls,vi_c1[,1], vi_h1[,1], vi_h2[,1]))
vi.melted <- melt(vi_data, id = "inter_ls")
vi.ci <- as.data.frame(cbind(inter_ls, vi_c1[,1] - (1.96*vi_c1[,2]/10), vi_c1[,1] + (1.96*vi_c1[,2]/10)))

#plot
vi_plot <- ggplot(data = vi.melted, aes(x = inter_ls, y = value,  color = variable)) +
  geom_point() + stat_smooth(se = F) + theme_classic() + 
  xlab("Interaction Strength") + ylab("% variation") + labs(title = "% variation in the outcome explained by interactions") + 
  geom_ribbon(data=vi.ci,aes(x=inter_ls,ymin=V2,ymax=V3),fill="slateblue2",alpha=0.2, inherit.aes = FALSE) + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)), legend.position="bottom") +
  scale_color_manual(values=c("sienna2","darkcyan", "midnightblue"), labels = c("1", "4", "10"), name="Heterogeneity Level")

show(vi_plot)


#(b) Repeated dataset
#First, run repeated.R

#vector of heterogeneity 'bad' level coefficients
het_ls <- seq(0, 10, .5)

#reading in files
means_c1s4 <- read.csv("means_c1s4")
sds_c1s4 <- read.csv("sds_c1s4")
means_c2s4 <- read.csv("means_c2s4")
sds_c2s4 <- read.csv("sds_c2s4")

#Creating the plotting data frames
cont_tcga <- as.matrix(cbind(means_c1s4$Unweighted,means_c1s4$Stack_ridge, means_c2s4$Stack_ridge))
data_tcga <- data.frame(cbind(het_ls, cont_tcga))
df.melted4 <- melt(data_tcga, id = "het_ls")
se_tcga.Tree <- data.frame(cbind(het_ls, means_c1s4$Stack_ridge - 1.96*sds_c1s4$Stack_ridge/(sqrt(100)),
                                 means_c1s4$Stack_ridge + 1.96*sds_c1s4$Stack_ridge/(sqrt(100))))
se_tcga.Forest <- data.frame(cbind(het_ls, means_c2s4$Stack_ridge - 1.96*sds_c2s4$Stack_ridge/(sqrt(100)),
                                   means_c2s4$Stack_ridge + 1.96*sds_c2s4$Stack_ridge/(sqrt(100))))
se_tcga.Unweighted <- data.frame(cbind(het_ls, means_c1s4$Unweighted - 1.96*sds_c1s4$Unweighted/(sqrt(100)),
                                       means_c1s4$Unweighted + 1.96*sds_c1s4$Unweighted/(sqrt(100))))

#Plot
q4 <- ggplot(data = df.melted4, aes(x = het_ls, y = value, color = variable)) +
  geom_point() + stat_smooth(se = F) + 
  geom_ribbon(data=se_tcga.Tree,aes(x=het_ls,ymin=V2,ymax=V3),fill="darkcyan",alpha=0.2, inherit.aes = FALSE) + 
  geom_ribbon(data=se_tcga.Forest,aes(x=het_ls,ymin=V2,ymax=V3),fill="midnightblue",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_tcga.Unweighted,aes(x=het_ls,ymin=V2,ymax=V3),fill="sienna2",alpha=0.2, inherit.aes = FALSE) +
  theme_classic() + xlab("Heterogeneity Level") + ylab("Percent change in average RMSE from Merged") +
  theme(legend.position="none") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5))) +
  labs(title = "Repeated dataset") +
  scale_color_manual(values=c("sienna2","darkcyan","midnightblue"), labels = c("Unweighted", "Weighting Trees", "Weighting Forests"), name="Method")
show(q4)


#(c) Splitting up the TCGA dataset
#First, run tcga_combined.R
het_ls <- seq(0, 10, .5)

#reading in files
means_c1s4 <- read.csv("means_c1s4")
sds_c1s4 <- read.csv("sds_c1s4")
means_c2s4 <- read.csv("means_c2s4")
sds_c2s4 <- read.csv("sds_c2s4")


#Creating the plotting data frames
cont_tcga <- as.matrix(cbind(means_c1s4$Unweighted,means_c1s4$Stack_ridge, means_c2s4$Stack_ridge))
data_tcga <- data.frame(cbind(het_ls, cont_tcga))
df.melted4 <- melt(data_tcga, id = "het_ls")
se_tcga.Tree <- data.frame(cbind(het_ls, means_c1s4$Stack_ridge - 1.96*sds_c1s4$Stack_ridge/(sqrt(100)),
                                 means_c1s4$Stack_ridge + 1.96*sds_c1s4$Stack_ridge/(sqrt(100))))
se_tcga.Forest <- data.frame(cbind(het_ls, means_c2s4$Stack_ridge - 1.96*sds_c2s4$Stack_ridge/(sqrt(100)),
                                   means_c2s4$Stack_ridge + 1.96*sds_c2s4$Stack_ridge/(sqrt(100))))
se_tcga.Unweighted <- data.frame(cbind(het_ls, means_c1s4$Unweighted - 1.96*sds_c1s4$Unweighted/(sqrt(100)),
                                       means_c1s4$Unweighted + 1.96*sds_c1s4$Unweighted/(sqrt(100))))


#Plot
q5 <- ggplot(data = df.melted4, aes(x = het_ls, y = value, color = variable)) +
  geom_point() + stat_smooth(se = F) + 
  geom_ribbon(data=se_tcga.Tree,aes(x=het_ls,ymin=V2,ymax=V3),fill="darkcyan",alpha=0.2, inherit.aes = FALSE) + 
  geom_ribbon(data=se_tcga.Forest,aes(x=het_ls,ymin=V2,ymax=V3),fill="midnightblue",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se_tcga.Unweighted,aes(x=het_ls,ymin=V2,ymax=V3),fill="sienna2",alpha=0.2, inherit.aes = FALSE) +
  theme_classic() + xlab("Heterogeneity Level") + ylab("Percent change in average RMSE from Merged") +
  labs(title = "TCGA study") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5))) +
  scale_color_manual(values=c("sienna2","darkcyan", "midnightblue"), labels = c("Unweighted", "Weighting Trees", "Weighting Forests"), name="Method")
show(q5)



#Make Figure 1
library(ggpubr)
ga <- ggarrange(q4, q5, labels = c("B", "C"), ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
vi_plot <- ggarrange(vi_plot, labels = c("A"), ncol = 1, nrow = 1)
ggarrange(vi_plot, ga, widths = c(1.1, 2), heights = c(1.7, 1))


