
#For A.1 - A.2, run multistudy_example_binary.R

#(A.1)
#Creating plotting data frame
plot_means <- data.frame(Method = c("Merged", "Unweighted", "Weighting Trees", "Weighting Forests"), 
                         Means = cbind(rep_test_tree$means[c(2, 13)], 
                                       rep_test_forest$means[c(1, 13)])[c(3,1,2,4)],
                         Sds = cbind(rep_test_tree$sds[c(2, 13)], 
                                     rep_test_forest$sds[c(1, 13)])[c(3,1,2,4)])

#plot
p1 <- ggplot(plot_means, aes(x=Method, y = Means, fill = Method)) + geom_bar(stat = "identity", fill = c("#E76E51","#F4A261","#254653", "#2A9D8F")) +
  theme_classic() + ylab("Average Log Loss") + labs(title = "Binary outcome (OS)") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)))

show(p1)

#(A.2) Zoomed in view, without the Merged 
#Creating plotting data frame
plot_means.z <- data.frame(Method = c("Unweighted", "Weighting Trees", "Weighting Forests"), 
                           Means = cbind(rep_test_tree$means[c(2, 13)], 
                                         rep_test_forest$means[c(1, 13)])[c(3,1,2,4)][-1],
                           Sds = cbind(rep_test_tree$sds[c(2, 13)], 
                                       rep_test_forest$sds[c(1, 13)])[c(3,1,2,4)][-1])

#plot
p2 <- ggplot(plot_means.z, aes(x=Method, y = Means, fill = Method)) + geom_bar(stat = "identity", fill = c("#F4A261","#254653", "#2A9D8F")) +
  theme_classic() + geom_errorbar(aes(ymin=Means-Sds/10, ymax=Means+Sds/10), width=.2) + 
  ylab("Average Log Loss") + labs(title = "Binary outcome (OS), zoomed in view") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)))

show(p2)



#For B, run multistudy_example_continuous.R
#(B)

#Reading in data
allmeans_tree <- read.csv("~/Desktop/real_data_continuous/allmeans_tree.csv", row.names=NULL)
allmeans_forest <- read.csv("~/Desktop/real_data_continuous/allmeans_forest.csv", row.names=NULL)


#Creating plotting dataframe
id <- rep(1:1000, 3)
allmeans <- as.data.frame(cbind(factor(id), as.factor(rep(c(1, 2, 3), each = 1000)), 
                                c(allmeans_tree$Unweighted, allmeans_tree$Stack_ridge, allmeans_forest$Stack_ridge)))

plot_means.c <- data.frame(Method = c("Unweighted", "Weighting Trees", "Weighting Forests"), 
                           Means = c(mean(allmeans_tree$Unweighted), mean(allmeans_tree$Stack_ridge), mean(allmeans_forest$Stack_ridge)),
                           Sds = c(sd(allmeans_tree$Unweighted), sd(allmeans_tree$Stack_ridge), 
                                   sd(allmeans_forest$Stack_ridge)))


#plot
c1 <- ggplot(plot_means.c, aes(x=Method, y = Means, fill = Method)) + geom_bar(stat = "identity", fill = c("#F4A261","#254653", "#2A9D8F")) +
  theme_classic() + geom_errorbar(aes(ymin=Means-Sds/100, ymax=Means+Sds/100), width=.2) + 
  ylab("Average Log Loss") + labs(title = "Binary outcome (OS), zoomed in view") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)))

show(c1)

#Creating Figure 5
library(gridExtra)
library(ggpubr)
ggarrange(p1, p2, c1, labels = c("A.1", "A.2", "B"), ncol=3, nrow=1)


