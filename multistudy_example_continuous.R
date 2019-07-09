install.packages("randomForest", repos="http://cran.r-project.org", dependencies = TRUE)
#install.packages("MLmetrics", repos="http://cran.r-project.org", dependencies = TRUE)
install.packages("glmnet", repos="http://cran.r-project.org", dependencies = TRUE)
library(MLmetrics)
library(randomForest)
#library(curatedOvarianData)
library(glmnet)
library(rpart)
library(nnet)
library(nnls)
library(foreach)
library(parallel)
library(doParallel)
library(genefilter)
library(abind)
library(BiocManager)



#source("http://bioconductor.org/biocLite.R")
#BiocManager::install("curatedBreastData")
library("curatedBreastData")
data(curatedBreastDataExprSetList)



#process only the first two datasets to avoid a long-running example:
#take top 5000 genes by variance from each dataset.
proc_curatedBreastDataExprSetList <- processExpressionSetList(exprSetList=curatedBreastDataExprSetList[10:15], outputFileDirectory = "./", numTopVarGenes=5000)

esets <- lapply(1:length(proc_curatedBreastDataExprSetList), function(x) exprs(proc_curatedBreastDataExprSetList[[x]]))


cn <- lapply(esets, rownames)
cn_int <- Reduce(intersect, cn)

edat <- vector("list", length(esets))

for(i in 1:length(esets)){
  edat[[i]] <- esets[[i]][cn_int,]
}


for(i in 1:length(edat)){
  edat[[i]] <- t(edat[[i]])
  edat[[i]] <- as.data.frame(edat[[i]])
}


edat_orig <- edat


################################


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

randomforestfit <- function(data, complexity, numtree, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = numtree, maxnodes = complexity, ...)
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

multistudysim_tot <- function(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight, index){
  
  ndat <- length(edat_orig)
  ntrain <- 4
  nvar <- 100
  modfit <- modfit
  modpred <- modpred
  
  edat <- edat_orig
  
  #for continuous example, leave this uncommented
  for(i in 1:length(edat)){
    nvar <- ncol(edat[[i]]) - 1
    y <- edat[[i]][,index]
    edat[[i]] <- cbind(y, edat[[i]][,-index])
    colnames(edat[[i]]) <- c("y", paste0("V", 1:nvar))
  }
  #########
  
  # Learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain]) #binds together the first ntrain rows into data frame 
  
  
  if (treeweight == TRUE){
    mod0 <- modfit(matstack, complexity, numtree*ntrain) #applies the training function
  }
  else{
    mod0 <- modfit(matstack, complexity, numtree)
  }
  mods <- allpreds <- vector("list", ntrain)
  if (treeweight == TRUE){
    mses <- matrix(NA, ntrain*numtree, ntrain)
  }
  else{
    mses <- matrix(NA, ntrain, ntrain)
  }
  
  if (treeweight == TRUE){
    for(i in 1: ntrain){
      mods[[i]] <- modfit(edat[[i]], complexity, numtree)
      preds <- lapply(edat[1:ntrain], function(x){
        modpred(mods[[i]], newdata = x[, -1], treeweight) 
      })
      new_pred <- vector("list", numtree)
      for (j in 1:numtree){
        new_pred[[j]] <- vector("list", ntrain)
      }
      for (k in 1:ntrain){
        preds_k <- preds[[k]]
        for (h in 1:numtree){
          new_pred[[h]][[k]] <- as.vector(preds_k[,h])
        }
      }
      
      for (j in 1:numtree){
        if (i == 1){k <- j}
        else{k <- (i - 1)*numtree + j}
        mses[k,] <- unlist(lapply(1:ntrain, function(i){
          x = edat[[i]]
          if (bin == TRUE){
            LogLoss(as.numeric(as.character(x[,"y"])), as.numeric(new_pred[[j]][[i]]))}
          else{mean((new_pred[[j]][[i]] - x[,"y"])^2)}
        }))
        curpreds <- lapply(new_pred[[j]], as.numeric)
        allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
      }}
    for (i in 1:ntrain){
      k <- (i-1)*numtree + 1 
      mses[,i][k:(k+numtree-1)] <- NA 
    }}
  else{
    for(i in 1:ntrain){
      mods[[i]] <- modfit(edat[[i]], complexity, numtree)
      preds <- lapply(edat[1:ntrain], function(x){
        modpred(mods[[i]], newdata = x[, -1], treeweight) 
      })
      mses[i,] <- unlist(lapply(edat[1:ntrain], function(x){#cross validation within the training set
        newdata = x[, -1]
        preds <- modpred(mods[[i]], newdata = newdata, treeweight) 
        if (bin == TRUE){
          LogLoss(as.numeric(as.character(x[,"y"])), as.numeric(preds))
        }
        else{mean((preds - x[,"y"])^2)}}
      ))
      
      curpreds <- lapply(preds, as.numeric)
      allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)}
    diag(mses) <- NA
  }
  
  # CS Weights
  tt <- apply(mses, 1, mean, na.rm = T) #removing the diagonal elements, takes the mean of each row 
  weights <- absnorm(sqrt(tt), max.norm = TRUE)
  nk <- unlist(lapply(edat, nrow)) #vector of number of rows in each dataset
  nwts <- absnorm(nk[1:ntrain])
  if (treeweight == TRUE){
    nwts <- rep(nwts, each = numtree)
  }
  # Regression: stacked (intercept and no intercept)
  predstack <- do.call(rbind, allpreds)
  coefs_stack_noint <- nnls::nnls(predstack, as.numeric(as.character(matstack$y)))$x
  coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), as.numeric(as.character(matstack$y)))$x
  coefs_stack_lasso <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 1, lower.limits = 0, intercept = T)))
  coefs_stack_ridge <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 0, lower.limits = 0, intercept = T)))
  #Just a safeguard against full collinearity, although I think we are OK with nnls now
  coefs_stack_noint[which(is.na(coefs_stack_noint))] <- 0
  coefs_stack_int[which(is.na(coefs_stack_int))] <- 0
  coefs_stack_lasso[which(is.na(coefs_stack_lasso))] <- 0
  coefs_stack_ridge[which(is.na(coefs_stack_ridge))] <- 0
  
  coefs_stack_noint_norm <- absnorm(coefs_stack_noint)
  
  # Regression: study-specific (intercept and no intercept)
  coefs_ss_noint <- mapply(function(x,y){nnls::nnls(y,as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_noint <- colMeans(do.call(rbind, coefs_ss_noint), na.rm = T)
  coefs_ss_int <- mapply(function(x,y){nnls::nnls(cbind(rep(1,nrow(y)),y),as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_int <- colMeans(do.call(rbind, coefs_ss_int), na.rm = T)
  coefs_ss_lasso <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 1, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_ridge <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 0, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_lasso <- colMeans(do.call(rbind, coefs_ss_lasso), na.rm = T)
  coefs_ss_ridge <- colMeans(do.call(rbind, coefs_ss_ridge), na.rm = T)
  coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
  coefs_ss_int[which(is.na(coefs_ss_int))] <- 0
  coefs_ss_lasso[which(is.na(coefs_ss_lasso))] <- 0
  
  coefs_ss_noint_norm <- absnorm(coefs_ss_noint)
  
  coefs_of_interest <- cbind(coefs_stack_int, coefs_stack_lasso, coefs_stack_ridge, coefs_ss_int, coefs_ss_lasso, coefs_ss_ridge)
  
  #how many coefficients are zeroed out by stacking vs. lasso, as well as how many of the zeroed out coefficients are shared between the two approaches
  zeros_summary <- t(matrix(NA, 6, 1))
  colnames(zeros_summary) <- c("stack_0", "stacklasso_0","SS_0", "SSlasso_0", "stack_intersect", "SS_intersect")
  stack_i0 <- sum(coefs_stack_int == 0)
  stack_l0 <- sum(coefs_stack_lasso == 0)
  ss_i0 <- sum(coefs_ss_int == 0)
  ss_l0 <- sum(coefs_ss_lasso == 0)
  stack_inter <- length(intersect(which(coefs_stack_int == 0), which(coefs_stack_lasso == 0)))
  ss_inter <- length(intersect(which(coefs_ss_int == 0), which(coefs_ss_lasso == 0)))
  zeros_summary[1,] <- c(stack_i0, stack_l0, ss_i0, ss_l0, stack_inter, ss_inter)
  
  outmat <- matrix(NA, ndat - ntrain, 14)
  colnames(outmat) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                        "Stack_noint", "Stack_noint_norm", "Stack_int",
                        "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  
  for(i in (ntrain + 1):ndat){
    merged <- modpred(mod0, newdata = edat[[i]][,-1], treeweight)	
    if (treeweight == TRUE){
      merged <- apply(merged, 2, as.numeric)
      allmod <- t(do.call(cbind,lapply(mods, modpred, newdata = edat[[i]][,-1], treeweight)))
      allmod <- apply(allmod, 2, as.numeric)
    }
    else{
      merged <- as.vector(sapply(merged, as.numeric))
      allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1], treeweight))
      allmod <- apply(allmod, 2, as.numeric)
    }
    # Unweighted average
    unweighted <- colMeans(allmod)
    
    # sample size weighted
    sample_wtd <- apply(allmod, 2, function(x){sum(nwts*x)})
    
    # cross-study weighted 
    cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})
    
    # regression: stacked (noint, int, each normed) + lasso
    stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
    stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
    stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})
    stack_lasso <- apply(allmod, 2, function(x){coefs_stack_lasso[1] + sum(coefs_stack_lasso[-1]*x)})
    stack_ridge <- apply(allmod, 2, function(x){coefs_stack_ridge[1] + sum(coefs_stack_ridge[-1]*x)})
    
    # regression: study_specific (noint, int, noint normed) + lasso
    ss_noint <- apply(allmod, 2, function(x){sum(coefs_ss_noint*x)})
    ss_noint_norm <- apply(allmod, 2, function(x){sum(coefs_ss_noint_norm*x)})
    ss_int <- apply(allmod, 2, function(x){coefs_ss_int[1] + sum(coefs_ss_int[-1]*x)})
    ss_lasso <- apply(allmod, 2, function(x){coefs_ss_lasso[1] + sum(coefs_ss_lasso[-1]*x)})
    ss_ridge <- apply(allmod, 2, function(x){coefs_ss_ridge[1] + sum(coefs_ss_ridge[-1]*x)})
    
    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    
    if (bin == TRUE){
      outmat[i - ntrain,] <- c(LogLoss(merged, cury), LogLoss(unweighted, cury),  LogLoss(sample_wtd, cury), LogLoss(cs_wtd, cury),
                               LogLoss(stack_noint, cury), LogLoss(stack_noint_norm, cury), LogLoss(stack_int, cury), 
                               LogLoss(ss_noint, cury), LogLoss(ss_noint_norm, cury),
                               LogLoss(ss_int, cury), LogLoss(stack_lasso, cury), LogLoss(ss_lasso, cury), 
                               LogLoss(stack_ridge, cury), LogLoss(ss_ridge, cury))}
    else{
      outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - sample_wtd)^2), mean((cury - cs_wtd)^2),
                                    mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
                                    mean((cury - ss_noint)^2), mean((cury - ss_noint_norm)^2), 
                                    mean((cury - ss_int)^2), mean((cury - stack_lasso)^2), mean((cury - ss_lasso)^2), 
                                    mean((cury - stack_ridge)^2), mean((cury - ss_ridge)^2)))}
  }
  outmat <- (outmat - outmat[,1])/outmat[,1]*100
  ##### Uncomment this section if want to output individual testing point errors: for continuous outcome
  #print(colMeans(outmat))
  #return(list(true_outcome = cury, merged = rowMeans(cury - merged), unweighted = (cury - unweighted), stack = (cury - stack_int), ss = (cury - ss_int), stack_lasso = (cury - stack_lasso), ss_lasso = (cury - ss_lasso)))
  return(list(means = colMeans(outmat), zeros_summary = zeros_summary, coefs_list = coefs_of_interest))
}

#Tests 
#Weighting trees
#multistudysim_tot(modfit = randomforestfit, modpred = randomforestpredict, good = .25, bad = 1, val = .4, edat_orig, simtype = "normal", complexity = 10, numtree = 10, bin = FALSE, ninter = 0, icoefs = c(4.4, -1.8), setnc = TRUE, treeweight = TRUE, index = 1000)

#Weighting forests
#multistudysim_tot(modfit = randomforestfit, modpred = randomforestpredict, good = .25, bad = 1, val = .4, edat_orig, simtype = "normal", complexity = 10, numtree = 10, bin = FALSE, ninter = 0, icoefs = c(4.4, -1.8), setnc = TRUE, treeweight = FALSE, index = 5)

 
#Now, replicate this 100 times each
rep_multistudy <- function(max_index, modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight, allmeansfile, coeffile){
  ptm <- proc.time()
  total_errors <- matrix(NA, max_index, 14)
  zeros_summary <- matrix(NA, max_index, 6)
  colnames_totalerror <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                           "Stack_noint", "Stack_noint_norm", "Stack_int",
                           "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  colnames(total_errors) <- colnames_totalerror
  colnames(zeros_summary) <- c("stack_0", "stacklasso_0", "SS_0", "SSlasso_0", "stack_intersect", "SS_intersect")
  
  for (i in 1:max_index){
    mssi <- c(multistudysim_tot(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight, i))
    total_errors[i, ] <- mssi$means
    zeros_summary[i, ] <- mssi$zeros_summary
    if (i == 1){
      coefs_out <- mssi$coefs_list
    }
    else{
      coefs_out <- rbind(coefs_out, mssi$coefs_list)
    }
    
    if (i %% 10 == 0){
      cat(paste0("iteration: ", i, "\n"))
    }
  }
  coefs_out <- as.data.frame(coefs_out)
  coef_colnames <- c("coefs_stack_int", "coefs_stack_lasso", "coefs_stack_ridge", "coefs_ss_int", "coefs_ss_lasso", "coefs_ss_ridge")
  colnames(coefs_out) <- coef_colnames
  total_errors <- na.omit(total_errors)
  zeros_summary <- colMeans(zeros_summary)
  write.table(total_errors, allmeansfile, sep = ",", row.names = F, col.names = colnames_totalerror)
  write.table(coefs_out, coeffile, sep = ",", row.names = F, col.names = coef_colnames)
  means <- colMeans(total_errors)
  sds <- apply(total_errors, 2, sd)
  write.table(t(as.matrix(means)), paste0("means", treeweight), sep = ",", row.names = F, col.names = colnames_totalerror)
  write.table(t(as.matrix(sds)), paste0("sds", treeweight), sep = ",", row.names = F, col.names = colnames_totalerror)
  print(proc.time() - ptm)
  return(list(total_errors = total_errors, zeros = zeros_summary, total_coefs = coefs_out, means = means, sds = sds)) 
}


#Weighting forests
rep_test_forest <- rep_multistudy(max_index = 1000, modfit = randomforestfit, modpred = randomforestpredict, good = .25, bad = 1, val = .4, edat_orig, simtype = "normal", complexity = 10, numtree = 10, bin = FALSE, ninter = 0, icoefs = c(4.4, -1.8), setnc = TRUE, treeweight = FALSE, allmeansfile = "allmeans_forest.csv", coeffile = "coefs_forest.csv")

#Weighting trees
rep_test_tree <- rep_multistudy(max_index = 1000, modfit = randomforestfit, modpred = randomforestpredict, good = .25, bad = 1, val = .4, edat_orig, simtype = "normal", complexity = 10, numtree = 10, bin = FALSE, ninter = 0, icoefs = c(4.4, -1.8), setnc = TRUE, treeweight = TRUE, allmeansfile = "allmeans_tree.csv", coeffile = "coefs_tree.csv")

##Plot
allmeans_tree <- read.csv("~/Desktop/real_data_continuous/allmeans_tree.csv", row.names=NULL)
allmeans_forest <- read.csv("~/Desktop/real_data_continuous/allmeans_forest.csv", row.names=NULL)

id <- rep(1:1000, 3)
allmeans <- as.data.frame(cbind(factor(id), as.factor(rep(c(1, 2, 3), each = 1000)), 
                                c(allmeans_tree$Unweighted, allmeans_tree$Stack_ridge, allmeans_forest$Stack_ridge)))

colnames(allmeans) <- c("Unweighted", "Weighting Trees", "Weighting Forests")


meansTrue <- read.csv("~/Desktop/real_data_continuous/meansTrue", row.names=NULL)
meansFalse <- read.csv("~/Desktop/real_data_continuous/meansFalse", row.names=NULL)
sdsTrue <- read.csv("~/Desktop/real_data_continuous/sdsTrue", row.names=NULL)
sdsFalse <- read.csv("~/Desktop/real_data_continuous/sdsFalse", row.names=NULL)


plot_means.c <- data.frame(Method = c("Unweighted", "Weighting Trees", "Weighting Forests"), 
                           Means = c(1.979604, 1.821666, 1.850352),
                           Sds = c(sdsTrue[c(2, 13)], 
                                       sdsFalse[c(1, 13)])[c(3,1,2,4)][-1])


plot_means.c <- data.frame(Method = c("Unweighted", "Weighting Trees", "Weighting Forests"), 
                         Means = c(mean(allmeans_tree$Unweighted), mean(allmeans_tree$Stack_ridge), mean(allmeans_forest$Stack_ridge)),
                         Sds = c(sd(allmeans_tree$Unweighted), sd(allmeans_tree$Stack_ridge), 
                                     sd(allmeans_forest$Stack_ridge)))

c1 <- ggplot(plot_means.c, aes(x=Method, y = Means, fill = Method)) + geom_bar(stat = "identity", fill = c("#F4A261","#254653", "#2A9D8F")) +
  theme_classic() + geom_errorbar(aes(ymin=Means-Sds/100, ymax=Means+Sds/100), width=.2) + 
  ylab("Average Log Loss") + labs(title = "Binary outcome (OS), zoomed in view") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5)))



colnames(allmeans) <- c("id", "method", "difference")

c1 <- ggplot(allmeans, aes(x=method, y=difference, fill = factor(method))) + 
  geom_boxplot(outlier.size = 1) + theme_classic() 

show(c1)
