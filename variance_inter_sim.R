library(MLmetrics)
library(randomForest)
library(curatedOvarianData)
library(glmnet)
library(rpart)
library(nnet)
library(nnls)
library(foreach)
library(parallel)
library(doParallel)
library(genefilter)
library(abind)
library(caret)
library(Biobase)
library(dplyr)

#Requires to have first saved the ovarian cancer datasets: in this case, it is saved as edat_orig.Rda
#For home computer
#load("~/Desktop/edat_orig.Rda")

#For cluster
load("/home/mr356/edat_orig.Rda")

# Work with the intersection of rows
cn <- lapply(esets, rownames)
cn_int <- Reduce(intersect, cn)

edat <- vector("list", length(esets))

for(i in 1:length(esets)){
  edat[[i]] <- esets[[i]][cn_int,]
}

# Remove esets with missing gene expression data
ridx <- which(unlist(lapply(edat, function(x){sum(is.na(exprs(x)))})) > 0)
edat <- edat[-ridx] # length 15

# Convert eset list to set of matrices
for(i in 1:length(edat)){
  edat[[i]] <- t(exprs(edat[[i]]))
}

edat_orig <- edat

# Normalize the columns
for(i in 1:length(edat_orig)){
  edat_orig[[i]] <- apply(edat_orig[[i]], 2, scale)
}



randomforestfit <- function(data, complexity, numtree, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = numtree, maxnodes = complexity, importance = TRUE, ...)
  rf
}

randomforestpredict <- function(data, newdata, treeweight){
  #newdata = as.data.frame(newdata)
  if (treeweight == TRUE){
    pred_obj <- predict(data, newdata=newdata, predict.all = TRUE)
    pred_obj$individual
  }
  else{
    as.vector(predict(data, newdata=newdata))
  }
}


#Error function for binary data
LogLoss <- function (y_pred, y_true){
  eps <- 1e-15
  y_pred <- pmax(pmin(y_pred, 1 - eps), eps)
  LogLoss <- -mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred))
  return(LogLoss)
}

absnorm <- function(vec, max.norm = FALSE){
  sgn <- sign(vec)
  vec <- abs(vec)
  if(max.norm){
    mvec <- max(vec)	
  } else {
    mvec <- min(vec)
  }
  
  vec <- abs(vec - mvec)
  sgn*(vec/sum(vec))
}



#update in order to output which variables are used to generate the outcome
init_data <- function(edat_orig, ndat, nvar, simtype, ninter, good, bad, val, icoefs, setnc, bin){
  
  edat <- edat_orig
  edat <- edat[sample(1:ndat)] 
  
  idx <- sample(1:ncol(edat[[1]]), nvar)
  for(i in 1:ndat){
    edat[[i]] <- edat[[i]][,idx] 
  }
  
  if (setnc == TRUE){ #setting the number of coefficients used to generate outcome
    ncoef <- 10} #hardcoded: need to change this number
  else{#random number of coefficients used to generate the outcome
    #Generate Y - outcome for each patient given the genes selected 
    if(simtype == "nonl"){
      ncoef <- sample(3:nvar, 1)
    } else {
      ncoef <- sample(2:nvar, 1) #number of coefficients used to create the outcome is between 2 and nvar, sampled randomly
    }}
  
  coefs <- sample(c(runif(round(ncoef/2), -5, -0.5), runif(ncoef - round(ncoef/2), 0.5, 5))) 
  vars <- sample(1:ncol(edat[[1]]), ncoef) 
  
  #the mean percentage of the outcome (y) that is explained by the interaction terms
  #num_columns = number of interactions
  variance_inter <- 0
  for(i in 1:ndat){ #ndat is between 1 and 15 (number of datasets)
    
    curcoefs <- sapply(coefs, function(x){runif(1, x - (i<=5)*good - 
                                                  (i > 5 & i <= 10)*bad - (i > 10)*val , x + 
                                                  (i<=5)*good + (i > 5 & i <= 10)*bad + (i > 10)*val)}) #adds noise to the coefficients 
    
    
    if(simtype == "slash"){
      y <- (edat[[i]][,vars] %*% curcoefs) + 
        0.05*cbind(rnorm(nrow(edat[[i]]))/runif(nrow(edat[[i]]))) # Added "slash" noise
    } else if(simtype == "nonl"){
      if ((1 <= i) & (i <= ninter)){
        y <- (edat[[i]][,vars] %*% curcoefs) + icoefs[1]*edat[[i]][,vars[1]]*edat[[i]][,vars[2]] 
        - icoefs[2]*edat[[i]][,vars[1]]*edat[[i]][,vars[3]] + cbind(rnorm(nrow(edat[[i]]))) # Added interaction terms
        #variance_inter_i <- (icoefs[1]*edat[[i]][,vars[1]]*edat[[i]][,vars[2]] 
        #              - icoefs[2]*edat[[i]][,vars[1]]*edat[[i]][,vars[3]])/y
        #variance_inter <- c(variance_inter, mean(variance_inter_i) * 100)
        xmat <- as.data.frame(edat[[i]][,vars])
        xmat <- cbind(y, xmat)
        #xmat_inter <- as.data.frame(cbind(edat[[i]][,vars], edat[[i]][,vars[1]]*edat[[i]][,vars[2]], edat[[i]][,vars[1]]*edat[[i]][,vars[3]]))
        #xmat_inter <- cbind(y, xmat_inter)
        lm_nointer <- lm(y ~ ., data = xmat)
        #lm_inter <- lm(y ~., data = xmat_inter)
        r2.lm_ni <- summary(lm_nointer)$r.squared
        #r2.lm_i <- summary(lm_inter)$r.squared
        #partial.r2 <- r2.lm_i- r2.lm_ni
        partial.r2 <- 1- r2.lm_ni
        variance_inter <- variance_inter + ((1/ninter)*partial.r2)
      }
      else{
        y <- (edat[[i]][,vars] %*% curcoefs) + cbind(rnorm(nrow(edat[[i]]))) # Added noise
      }
    } else {
      y <- (edat[[i]][,vars] %*% curcoefs) + cbind(rnorm(nrow(edat[[i]]))) # Added noise
    }
    
    
    if (bin == TRUE){
      new_y <- ifelse(y >= quantile(y, .75), 1, 0)
      
      edat[[i]] <- cbind(new_y, edat[[i]])
      edat[[i]] <- as.data.frame(edat[[i]])
      colnames(edat[[i]]) <- c("y", paste0("V", 1:nvar))
    }
    else{
      edat[[i]] <- cbind(y, edat[[i]])
      edat[[i]] <- as.data.frame(edat[[i]])
      colnames(edat[[i]]) <- c("y", paste0("V", 1:nvar))
    }
  }
  if (bin == TRUE){
    for(i in 1:ndat){
      edat[[i]]$y <- factor(edat[[i]]$y)
    }
  }
  return(list(edat = edat, idx = idx, vars = vars, curcoefs = curcoefs, variance_inter = variance_inter))
}

vary_init <- function(reps, var_list, edat_orig, ndat, nvar, simtype, ninter, good, bad, val, icoefs, setnc, bin){
  ptm = proc.time()
  total_variance_inter <- array(0, c(length(var_list),2))
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    icoefs <- c(level, -level/3)
    total_variance_inter_i <- array(0, c(reps, 2))
    for (j in 1:reps){
      id_test <- init_data(edat_orig, ndat, nvar, simtype, ninter, good, bad, val, icoefs, setnc, bin)
      total_variance_inter_i[j,] <- id_test$variance_inter
    }
    total_variance_inter[i,] <- c(mean(total_variance_inter_i), sd(total_variance_inter_i))
  }
  return(total_variance_inter)
}

inter_ls <- seq(2, 12, .5)

vi_c1 <- vary_init(reps = 100, var_list = inter_ls,edat_orig, ndat = length(edat_orig), nvar = 100, simtype = "nonl", ninter = 4, good = .25, bad = 1, val = .4, icoefs = c(4.4, -1.8), setnc = TRUE, bin = FALSE)
#vi_c2 <- vary_init(reps = 100, var_list = inter_ls,edat_orig, ndat = length(edat_orig), nvar = 100, simtype = "nonl", ninter = 8, good = .25, bad = 1, val = .4, icoefs = c(4.4, -1.8), setnc = TRUE, bin = FALSE)


#change heterogeneity levels
#h1: level = 4
vi_h1 <- vary_init(reps = 100, var_list = inter_ls, edat_orig, ndat = length(edat_orig), nvar = 100, simtype = "nonl", ninter = 4, good = .25, bad = 4, val = 2, icoefs = c(4.4, -1.8), setnc = TRUE, bin = FALSE)
#h2: level = 10
vi_h2 <- vary_init(reps = 100, var_list = inter_ls, edat_orig, ndat = length(edat_orig), nvar = 100, simtype = "nonl", ninter = 4, good = .25, bad = 10, val = 5, icoefs = c(4.4, -1.8), setnc = TRUE, bin = FALSE)


#plot
#vi_data <- as.data.frame(cbind(inter_ls, vi_c1[,1], vi_c1[,1] - (1.96*vi_c1[,2]/10), vi_c1[,1] + (1.96*vi_c1[,2]/10)))
vi_data <- as.data.frame(cbind(inter_ls,vi_c1[,1], vi_h1[,1], vi_h2[,1]))
vi.melted <- melt(vi_data, id = "inter_ls")
vi.ci <- as.data.frame(cbind(inter_ls, vi_c1[,1] - (1.96*vi_c1[,2]/10), vi_c1[,1] + (1.96*vi_c1[,2]/10)))

vi_plot <- ggplot(data = vi.melted, aes(x = inter_ls, y = value,  color = variable)) +
  geom_point() + stat_smooth(se = F) + theme_classic() + 
  xlab("Interaction Strength") + ylab("% variation") + labs(title = "% variation in the outcome explained by interactions") + 
  geom_ribbon(data=vi.ci,aes(x=inter_ls,ymin=V2,ymax=V3),fill="slateblue2",alpha=0.2, inherit.aes = FALSE) + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)), legend.position="bottom") +
  scale_color_manual(values=c("sienna2","darkcyan", "midnightblue"), labels = c("1", "4", "10"), name="Heterogeneity Level")
  #scale_color_manual(values=c("#00AEDB","#00AFAF","#00B084"), labels = c("1", "4", "10"), name="Heterogeneity Level")
  
show(vi_plot)
