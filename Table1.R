#First, run main_simulation.R

##################################
#Tables for Variable importance 

#Case 1
varimp_c1_forest <- read.csv("varimp_c1_forest.csv", row.names=NULL)
varimp_c1_merged <- read.csv("varimp_c1_merged.csv", row.names=NULL)
varimp_c1_tree <- read.csv("varimp_c1_tree.csv", row.names=NULL)


library(dplyr)
#Forest
#top 10
vi1.forest.top.mean <- as.data.frame(varimp_c1_forest[c( TRUE, rep(FALSE, 2) ), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.forest.top.sd <- as.data.frame(varimp_c1_forest[c(TRUE, rep(FALSE, 2)), ] %>% group_by(row.names) %>% summarise_all(funs(sd)))

#bottom 10
vi1.forest.bottom.mean <- as.data.frame(varimp_c1_forest[c( FALSE, TRUE, FALSE ), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.forest.bottom.sd <- as.data.frame(varimp_c1_forest[c( FALSE, TRUE, FALSE ), ] %>% group_by(row.names) %>% summarise_all(funs(sd)))

#total
vi1.forest.total.mean <- as.data.frame(varimp_c1_forest[c( FALSE, FALSE, TRUE ), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.forest.total.sd <- as.data.frame(varimp_c1_forest[c( FALSE, FALSE, TRUE ), ] %>% group_by(row.names) %>% summarise_all(funs(sd)))

vi1.forest.top.mean <- round(colMeans(vi1.forest.top.mean[order(as.numeric(vi1.forest.top.mean$row.names)),][, -1]), digits = 3)
vi1.forest.bottom.mean <- round(colMeans(vi1.forest.bottom.mean[order(as.numeric(vi1.forest.bottom.mean$row.names)),][, -1]), digits = 3)
vi1.forest.total.mean <- round(colMeans(vi1.forest.total.mean[order(as.numeric(vi1.forest.total.mean$row.names)),][, -1]), digits = 3)

vi1.forest.means <- rbind(vi1.forest.top.mean, vi1.forest.bottom.mean, vi1.forest.total.mean)
vi1.forest.means <- cbind(vi1.forest.means, round(vi1.forest.means[, 6]/vi1.forest.means[, 7], digits = 3))[,c(1, 2, 3, 4, 9)]

vi1.forest.top.sd <- round(colMeans(vi1.forest.top.sd[order(as.numeric(vi1.forest.top.sd$row.names)),][, -1]), digits = 3)
vi1.forest.bottom.sd <- round(colMeans(vi1.forest.bottom.sd[order(as.numeric(vi1.forest.bottom.sd$row.names)),][, -1]), digits = 3)
vi1.forest.total.sd <- round(colMeans(vi1.forest.total.sd[order(as.numeric(vi1.forest.total.sd$row.names)),][, -1]), digits = 3)

vi1.forest.sd <- rbind(vi1.forest.top.sd, vi1.forest.bottom.sd, vi1.forest.total.sd)
vi1.forest.sd <- cbind(vi1.forest.sd, round(vi1.forest.sd[, 6]/vi1.forest.sd[, 7], digits =3))[,c(1, 2, 3, 4, 9)]


msd <- paste("\\vtop{\\hbox{\\strut \\textbf{",vi1.forest.means,"}}\\hbox{(",round(vi1.forest.means - 1.96*vi1.forest.sd/sqrt(2100), digits = 3), ", ", 
             round(vi1.forest.means + 1.96*vi1.forest.sd/sqrt(2100), digits = 3),")","}}",sep="")

tab <- matrix(msd,3,5)[,-2]
#colnames <- c("Avg. tree weight", "# variables per tree", "freq. of true variables", "# variables in interactions", "sum of variable importance scores \n for true variables",
#              "sum of variable importance scores for all variables", "sum of true variable coefs")
#colnames(tab) <- colnames
options(xtable.sanitize.text.function=identity)
x = xtable(tab, digits = 3)
print(x, booktabs = getOption("xtable.booktabs", TRUE))





##################
#Tree
#top 10
vi1.tree.top.mean <- as.data.frame(varimp_c1_tree[c( TRUE, rep(FALSE, 2) ), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.tree.top.sd <- as.data.frame(varimp_c1_tree[c(TRUE, rep(FALSE, 2)), ] %>% group_by(row.names) %>% summarise_all(funs(sd)))

#bottom 10
vi1.tree.bottom.mean <- as.data.frame(varimp_c1_tree[c( FALSE, TRUE, FALSE ), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.tree.bottom.sd <- as.data.frame(varimp_c1_tree[c( FALSE, TRUE, FALSE ), ] %>% group_by(row.names) %>% summarise_all(funs(sd)))

#total
vi1.tree.total.mean <- as.data.frame(varimp_c1_tree[c( FALSE, FALSE, TRUE ), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.tree.total.sd <- as.data.frame(varimp_c1_tree[c( FALSE, FALSE, TRUE ), ] %>% group_by(row.names) %>% summarise_all(funs(sd)))

vi1.tree.top.mean <- round(colMeans(vi1.tree.top.mean[order(as.numeric(vi1.tree.top.mean$row.names)),][, -1]), digits = 3)
vi1.tree.bottom.mean <- round(colMeans(vi1.tree.bottom.mean[order(as.numeric(vi1.tree.bottom.mean$row.names)),][, -1]), digits = 3)
vi1.tree.total.mean <- round(colMeans(vi1.tree.total.mean[order(as.numeric(vi1.tree.total.mean$row.names)),][, -1]), digits = 3)

vi1.tree.means <- rbind(vi1.tree.top.mean, vi1.tree.bottom.mean, vi1.tree.total.mean)
vi1.tree.means <- cbind(vi1.tree.means, round(vi1.tree.means[, 6]/vi1.tree.means[, 7], digits = 3))[,c(1, 2, 3, 4, 9)]

vi1.tree.top.sd <- round(colMeans(vi1.tree.top.sd[order(as.numeric(vi1.tree.top.sd$row.names)),][, -1]), digits = 3)
vi1.tree.bottom.sd <- round(colMeans(vi1.tree.bottom.sd[order(as.numeric(vi1.tree.bottom.sd$row.names)),][, -1]), digits = 3)
vi1.tree.total.sd <- round(colMeans(vi1.tree.total.sd[order(as.numeric(vi1.tree.total.sd$row.names)),][, -1]), digits = 3)

vi1.tree.sd <- rbind(vi1.tree.top.sd, vi1.tree.bottom.sd, vi1.tree.total.sd)
vi1.tree.sd <- cbind(vi1.tree.sd, round(vi1.tree.sd[, 6]/vi1.tree.sd[, 7], digits =3))[,c(1, 2, 3, 4, 9)]

#msd <- paste(vi1.tree.means," (",round(vi1.tree.means - 1.96*vi1.tree.sd/sqrt(2100), digits = 3), ",", 
#                                   round(vi1.tree.means + 1.96*vi1.tree.sd/sqrt(2100), digits = 3),")",sep="")
msd <- paste("\\vtop{\\hbox{\\strut \\textbf{",vi1.tree.means,"}}\\hbox{(",round(vi1.tree.means - 1.96*vi1.tree.sd/sqrt(2100), digits = 4), ", ", 
             round(vi1.tree.means + 1.96*vi1.tree.sd/sqrt(2100), digits = 4),")","}}",sep="")

tab <- matrix(msd,3,5)[,-2]
#colnames <- c("Avg. tree weight", "# variables per tree", "freq. of true variables", "# variables in interactions", "sum of variable importance scores \n for true variables",
#              "sum of variable importance scores for all variables", "sum of true variable coefs")
#colnames(tab) <- colnames
options(xtable.sanitize.text.function=identity)
x = xtable(tab, digits = 4)
print(x, booktabs = getOption("xtable.booktabs", TRUE))


##################
#Merged 
varimp_c1_merged
#means
vi1.merged.group.mean <- as.data.frame(varimp_c1_merged[c( TRUE, FALSE ), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.merged.means <- round(colMeans(vi1.merged.group.mean[order(as.numeric(vi1.merged.group.mean$row.names)),][, -1]), digits = 3)
vi1.merged.means <- cbind(vi1.merged.means, round(vi1.merged.means[6]/vi1.merged.means[7], digits = 3))[c(1, 2, 3, 4, 9)]

#sds 
vi1.merged.group.sd <- as.data.frame(varimp_c1_merged[c(FALSE,TRUE), ] %>% group_by(row.names) %>% summarise_all(funs(mean)))
vi1.merged.sd <- round(colMeans(vi1.merged.group.sd[order(as.numeric(vi1.merged.group.sd$row.names)),][, -1]), digits = 3)
vi1.merged.sd <- cbind(vi1.merged.sd, round(vi1.merged.sd[6]/vi1.merged.sd[7], digits = 3))[c(1, 2, 3, 4, 9)]


#msd <- paste(vi1.merged.means," (",vi1.merged.means + 1.96*vi1.merged.sd/sqrt(2100),")",sep="")
msd <- paste("\\vtop{\\hbox{\\strut \\textbf{",vi1.merged.means,"}}\\hbox{(",round(vi1.merged.means - 1.96*vi1.merged.sd/sqrt(2100), digits = 3), ", ", 
             round(vi1.merged.means + 1.96*vi1.merged.sd/sqrt(2100), digits = 3),")","}}",sep="")

tab <- matrix(msd,1,5)[,-2]
#colnames <- c("Avg. tree weight", "# variables per tree", "freq. of true variables", "# variables in interactions", "sum of variable importance scores \n for true variables",
#              "sum of variable importance scores for all variables", "sum of true variable coefs")
#colnames(tab) <- colnames
options(xtable.sanitize.text.function=identity)
x = xtable(t(as.matrix(tab)), digits = 3)
print(x, booktabs = getOption("xtable.booktabs", TRUE))

