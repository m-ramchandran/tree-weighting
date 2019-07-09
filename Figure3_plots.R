#First, run main_simulation.R

#Read in the data
#Coefficients for level of het = 4, 2 interactions in training and 2 interactions in testing
coefs_c1_tree <- read.csv("coefs_c1s2_4.csv_tree.csv")
coefs_c1_forest <- read.csv("coefs_c1s2_4.csv_forest.csv")

#remove intercept terms
index_remove <- sapply(1:100, function(i) (i-1)*101 + 1)

#(a) Overall distribution of weights

#creating plotting data framces
coef_mat <- c(coefs_c1_forest[,3][-index_remove], coefs_c1_tree[,3][-index_remove], coefs_c1_forest[,3][-index_remove] - coefs_c1_tree[,3][-index_remove])
col_coef <- as.factor(rep(c(1, 2, 3), each = 10100 - length(index_remove)))
coef_mat <- as.data.frame(cbind(factor(col_coef), coef_mat))

#plot
i1 <- ggplot(coef_mat, aes(x=col_coef, y=coef_mat, fill = col_coef)) + 
  geom_boxplot(outlier.size = 1) + theme_classic()  + #+ labs(title = "Distribution of Weights")
  ylab("Tree weight") + scale_x_discrete(labels=c("1" = "Weighting Forests", "2" = "Weighting Trees", "3" = "Difference")) +
  scale_fill_brewer(palette="PuBuGn") + xlab("Method") + theme(legend.position="none") + 
  scale_y_continuous(breaks = seq(-.15, .15, .15), limits = c(-.15, .15)) +
  theme(axis.title=element_text(size=rel(1.7)), axis.text=element_text(size=rel(1.6)), plot.title = element_text(size = rel(1.5))) +
  labs(title = "Distribution of all weights") +
  #geom_hline(yintercept=.01, linetype="dashed", color = "red") + 
  geom_segment( aes(x=.3, y=0.01, xend=2.4,yend=0.01), linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0, linetype = "longdash", color = "black")
show(i1)

#(b) Distribution of weights for datasets with interaction terms 

#Indices of trees trained on datasets with interactions: 
index_with_interaction <- c()
for (j in 1:100){
  start <- (j-1)*100 + 1
  index_with_interaction <- c(index_with_interaction, start:(start + 19))
}

#create plotting data frames
coef_mat_inter <- c(coefs_c1_forest[,3][-index_remove][index_with_interaction], coefs_c1_tree[,3][-index_remove][index_with_interaction], coefs_c1_forest[,3][-index_remove][index_with_interaction] - coefs_c1_tree[,3][-index_remove][index_with_interaction])
col_coef_inter <- as.factor(rep(c(1, 2, 3), each = 2000))
coef_mat_inter <- as.data.frame(cbind(factor(col_coef_inter), coef_mat_inter))

i2 <- ggplot(coef_mat_inter, aes(x=col_coef_inter, y=coef_mat_inter, fill = col_coef_inter)) + 
  geom_boxplot(outlier.size = 1) + theme_classic()  + #+ labs(title = "Distribution of Weights")
  ylab("Tree weight") + scale_x_discrete(labels=c("1" = "Weighting Forests", "2" = "Weighting Trees", "3" = "Difference")) +
  scale_fill_brewer(palette="PuBuGn") + xlab("Method") + theme(legend.position="none") + 
  theme(axis.title=element_text(size=rel(1.7)), axis.text=element_text(size=rel(1.6)), plot.title = element_text(size = rel(1.5))) +
  labs(title = "Distribution of weights for trees trained on datasets with interactions") +  
  geom_hline(yintercept = 0, linetype = "longdash", color = "black") +
  geom_segment( aes(x=.3, y=0.01, xend=2.4,yend=0.01), linetype = "dashed", colour = "red") +
  #geom_hline(yintercept=.01, linetype="dashed", color = "red") + 
  scale_y_continuous(breaks = seq(-.15, .15, .15), limits = c(-.15, .15)) 

show(i2)


#Create Figure 3
library(gridExtra)
library(ggpubr)
ggarrange(i1,i2, labels = c("A", "B"), ncol=2, nrow=1)


