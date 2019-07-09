#Install curatedOvarianData
source("https://bioconductor.org/biocLite.R")
biocLite("curatedOvarianData")

# From curatedOvarianData vignette
source(system.file("extdata", "patientselection.config",package="curatedOvarianData"))
sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

# Save the original eset list - reuse this list for other analyses
# e.g. save(esets, file = "061417_esets.Rda")


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


###########################
#Helper functions for assessing variable importance and tree structure measures

#9 is the max number of variables per tree
#prints out the variables found in the trees of that forest
#each row represents the variables from that tree
#if m = mods:
#for the ith forest (i = 1,...10)
#v = 1:10
#if m = mod0:
#for the merged (i = 0)
#v = 1:100
var_modi <- function(i, v, m) {
  if (i == 0){
    ltest_i <- lapply(v, function(x) as.vector(na.omit(getTree(m, k = x, labelVar = TRUE)$`split var`)))
  }
  else{
    ltest_i <- lapply(v, function(x) as.vector(na.omit(getTree(m[[i]], k = x, labelVar = TRUE)$`split var`)))
  }
  ls <- do.call(rbind, lapply(v, function(y) c(ltest_i[[y]], rep(NA, 9 - length(ltest_i[[y]])))))
  row.names(ls) <- sapply(v, function(x) paste0("tree",x))
  colnames(ls) <- sapply(seq(1, ncol(ls)), function(x) paste0("var",x))
  return(ls)
}

#variable importance scores across all trees in forest i, if m = mods, or just in mod0, if m = mod0
#if i = 0, m = mod0
# if i = 1,...10, m = mods
sort_varimp_i <- function(i, m){
  if (i == 0){
    #vi <- varImp(m)$'0'
    vi <- varImp(m)
  }
  else{
    #vi <- varImp(m[[i]])$'0'
    vi <- varImp(m[[i]])
  }
  vi_order <- order(vi, decreasing = TRUE)
  varnames <- sapply(vi_order, function(x) paste0('V', x))
  newdata <- data.frame(matrix(ncol = 2, nrow = length(vi_order)))
  newdata[, 1] <- varnames
  #newdata[, 2] <- vi[vi_order]
  newdata[, 2] <- vi[vi_order,]
  colnames(newdata) <- c("Varname", "Importance")
  return(newdata)
}


#v = 1:10 
#m = mod0 or mods 
#m_ind = 0 if m = mod0, m_ind = 1 if m = mods
varImp_summary <- function(m, m_ind, vars_names, id, coefs_of_interest, simtype, v, numtree, ntrain){
  ntree = numtree*ntrain #total number of trees in ensemble
  if (m_ind == 1){
    #The variables contained in all 100 trees, when using SSL's (mods)
    mods_vars <- do.call(rbind, lapply(v, function(i) var_modi(i, v, m)))
    row.names(mods_vars) <- sapply(1:ntree, function(x) paste0("tree",x))
    
    #The sorted importance scores given to each variable for the 10 different forests (list where each index is the variable importance scores for that forest)
    mods_sort_varimp <- lapply(v, function(i) sort_varimp_i(i, m))
  }
  else{
    #The variables contained in all 100 trees, when using the merged (mod0)
    mods_vars <- var_modi(0, 1:ntree, m)
    
    #The sorted importance scores given to each variable for the 10 different forests (list where each index is the variable importance scores for that forest)
    mods_sort_varimp <- lapply(rep(0, ntrain), function(i) sort_varimp_i(0, m))
  }
  #list of lists; containing the variables in the intersection in that tree and the true variables
  tree_truth_mods <- lapply(1:ntree, function(x) intersect(mods_vars[x, ], vars_names)) 
  
  #vector; each index gives the proportion of that tree containing true variables
  freq_truth <- sapply(1:ntree, function(x) round(length(tree_truth_mods[[x]])/length(na.omit(mods_vars[x, ])), digits = 2)) 
  
  #vector; each index gives the number of variables in that tree
  num_var <- sapply(1:ntree, function(x) length(na.omit(mods_vars[x, ]))) 
  
  #The true coefficients given the the variables affecting the outcome 
  true_coefs_rank <- id$vars[order(id$curcoefs, decreasing = TRUE)]
  true_coefs_rank <- sapply(true_coefs_rank, function(i) paste0("V", i))
  true_coefs <- id$curcoefs[order(id$curcoefs, decreasing = TRUE)]
  names(true_coefs) <- true_coefs_rank
  
  #The true coefficients involved in interactions (if simtype == "nonl")
  inter_coefs <- true_coefs[1:3]
  
  #sum of the variable importance scores of the true variables contained in each tree
  sum_varimp <- c()
  
  #sum of the variable importance scores of all variables contained in each tree 
  sum_totalimp <- c()
  
  #Vector containing the sum of the coefficients (in the data generating mechanism) for the true variables contained in each tree
  sum_coef <- c()
  
  #Vector containing the number of true variables involved in interactions present
  num_inter <- c()
  
  #Vector containing the proportion of true variables included that are involved in interactions
  prop_inter <- c()
  
  for (i in 1:ntree){
    #To recover the forest that some tree (with index i) came from: 
    tree_index <- floor((i-1)/ntrain) + 1
    
    #a vector containing the importance scores for all variables in the forest that the i'th tree comes from:
    msvi <- mods_sort_varimp[[tree_index]]
    
    #True variable importance sum
    truevar_i <- tree_truth_mods[[i]] #the true variables contained in the i'th tree
    sum_varimp_i <- sum(msvi[match(truevar_i, msvi$Varname),]$Importance)
    sum_varimp <- c(sum_varimp, sum_varimp_i)
    
    #Total variable importance sum 
    total_var_i <- na.omit(mods_vars[i, ])
    sum_totalimp_i <- sum(msvi[match(total_var_i, msvi$Varname),]$Importance)
    sum_totalimp <- c(sum_totalimp, sum_totalimp_i)
    
    #sum of true coefs for each tree
    sum_coef_i <- sum(true_coefs[match(truevar_i, names(true_coefs))])
    sum_coef <- c(sum_coef, sum_coef_i)
    
    inter_indices_i <- match(names(inter_coefs), truevar_i)
    if ((length(inter_indices_i) == 1) && (is.na(inter_indices_i) == TRUE)){
      num_inter_i = 0 
      prop_inter_i = 0
    }
    else{
      inter_indices_i <- na.omit(inter_indices_i)
      num_inter_i <- length(inter_indices_i)
      if (length(truevar_i) == 0){
        prop_inter_i = 0
      }
      else{
        prop_inter_i <- num_inter_i/length(truevar_i) * 100
      }
    }
    num_inter <- c(num_inter, num_inter_i)
    prop_inter <- c(prop_inter, prop_inter_i)
  }
  
  
  if (m_ind == 1){
    #stack_ridge coefficients for every tree 
    #requires having run coefs of interest
    tree_coefs <- coefs_of_interest$coefs_stack_ridge[-1] #remove the intercept
  }
  else{
    tree_coefs <- rep(1/ntree, ntree)
  }
  #indices of trees, in order of weight
  order_stack <- order(tree_coefs, decreasing = TRUE)
  
  if (simtype == "nonl"){
    inter_coefs <- true_coefs[1:3]
    
    outmat <- cbind(tree_coefs, num_var, freq_truth, num_inter, prop_inter, sum_varimp, sum_totalimp, sum_coef)
  }
  else{
    outmat <- cbind(tree_coefs, num_var, freq_truth, sum_varimp, sum_totalimp, sum_coef)
  }
  
  #order in decreasing order of stack_ridge coefficients
  outmat <- as.data.frame(outmat[order_stack,])
  
  means_top10 <- round(apply(outmat[1:10,],2, mean), digits = 3)
  means_bottom10 <- round(apply(outmat[90:100,],2, mean), digits = 3)
  means_tot <- round(apply(outmat,2, mean), digits = 3)
  
  return(list(outmat = outmat, means_top10 = means_top10, means_bottom10 = means_bottom10, means_tot = means_tot))
}


#####################################
#fit and predict functions for Random Forest

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


# Set up a dataset split for one simulation iteration
#
# Input:
# edat_orig - orignal (cleaned) list of curatedOvarianData datasets
# ndat - length(edat_orig)
# nvar - number of features to reduce each dataset to
# simtype - one of "normal" (no interaction terms) or "nonl" (interaction terms)
# ninter - number of total datasets (both training and testing) with interaction terms
# good - Perturbation level for "good" i.e. low-perturbation studies
# bad - Perturbation level for "bad" i.e. high-perturbation studies
# val - Perturbation level for "val" i.e. validation studies
# icoefs - coefficients for interaction terms
# setnc - TRUE if set the number of coefficients used to generate outcome, FALSE if this number is randomly generated
# bin - TRUE if the outcome to be generated is binary, FALSE if it is continuous

# Output:
# edat - a list of esets with generated Y
# idx <- the indices of the variables included in each dataset from the original
# vars - the indices of the variables used to generate the outcome
# curcoefs - the values of the coefficients to generate the outcome from the chosen variables
# variance_inter - the amount of variance in the total outcome explained by interaction terms
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
  variance_inter <- c()
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
        variance_inter_i <- (icoefs[1]*edat[[i]][,vars[1]]*edat[[i]][,vars[2]] 
                             - icoefs[2]*edat[[i]][,vars[1]]*edat[[i]][,vars[3]])/y
        variance_inter <- c(variance_inter, mean(variance_inter_i) * 100)
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


############## Simulation functions

#Run once, then run multistudysim_tot with treeweight = TRUE and treeweight = FALSE, to ensure that the same trees are being weighted in each scheme and can be compared as such
multistudysim_helper <- function(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc){
  ndat <- length(edat_orig)
  ntrain <- 10
  nvar <- 100
  modfit <- modfit
  modpred <- modpred
  #this may need to change if the number of variables used to generate the outcome is not = 10
  v <- seq(1, numtree, by = 1)
  
  id <- init_data(edat_orig, ndat, nvar, simtype, ninter, good, bad, val, icoefs, setnc, bin)
  edat <- id$edat
  vars_names <- sapply(id$vars, function(x) paste0('V', x))
  variance_inter <- id$variance_inter
  
  if (ninter == 2){ #one study with interaction in training, one in testing
    edat <- replace(edat, c(2, 11), edat[c(11, 2)])
  } 
  else if (ninter == 4){
    numrepl <- 2
    edat <- replace(edat, (ninter - numrepl):ninter, edat[11:(11 + numrepl)])
  }
  else if (ninter > 4){
    numrepl <- floor(ninter/3)
    edat <- replace(edat, (ninter - numrepl):ninter, edat[11:(11 + numrepl)])
  }
  # Learn the merged predictor
  
  matstack <- do.call(rbind, edat[1:ntrain]) #binds together the first ntrain rows into data frame 
  
  mods <- vector("list", ntrain)
  for(i in 1: ntrain){
    mods[[i]] <- modfit(edat[[i]], complexity, numtree)
  }
  
  
  mod0 <- modfit(matstack, complexity, numtree*ntrain) #applies the training function
  
  
  #variable importance summary measures for mod0
  varimpsumm_mod0 <- varImp_summary(m = mod0, m_ind = 0, vars_names, id, coefs_of_interest, simtype, v, numtree, ntrain)
  
  return(list(id = id, edat = edat, vars_names = vars_names, variance_inter = variance_inter, matstack = matstack, mods = mods, mod0 = mod0, varimpsumm_mod0 = varimpsumm_mod0))
}


#Function for main simulation
#Input:
#ms_h - the result of running multistudysim_helper 
multistudysim_tot <- function(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight, ms_h){
  
  ndat <- length(edat_orig)
  ntrain <- 10
  nvar <- 100
  modfit <- modfit
  modpred <- modpred
  #this may need to change if the number of variables used to generate the outcome is not = 10
  v <- seq(1, numtree, by = 1)
  edat <- ms_h$edat
  vars_names <- ms_h$vars_names
  matstack <- ms_h$matstack
  mods <- ms_h$mods
  id <- ms_h$id
  mod0 <- ms_h$mod0
  variance_inter <- ms_h$variance_inter
  
  
  allpreds <- vector("list", ntrain)
  
  if (treeweight == TRUE){
    mses <- matrix(NA, ntrain*numtree, ntrain)
  }
  else{
    mses <- matrix(NA, ntrain, ntrain)
  }
  
  if (treeweight == TRUE){
    for(i in 1: ntrain){
      #mods[[i]] <- modfit(edat[[i]], complexity, numtree)
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
      #mods[[i]] <- modfit(edat[[i]], complexity, numtree)
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
  coefs_of_interest <- as.data.frame(coefs_of_interest)
  if (treeweight == FALSE){
    coefs_of_interest <- coefs_of_interest %>% slice(rep(1:n(), each = numtree))
    coefs_of_interest <- coefs_of_interest[-c(1:(numtree - 1)),]
    coefs_of_interest <- coefs_of_interest/10
    rownames(coefs_of_interest) <- 1:(ntrain*numtree + 1)
  }
  
  #Summary of variable importance measures for mods
  varimpsumm_mods <- varImp_summary(m = mods, m_ind = 1, vars_names, id, coefs_of_interest, simtype, v, numtree, ntrain)
  
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
  #Uncomment if want to return percent change from the merged:
  #outmat <- (outmat - outmat[,1])/outmat[,1]*100


  return(list(means = colMeans(outmat), coefs_list = coefs_of_interest, variance_inter = variance_inter, varimpsumm_mods = varimpsumm_mods))
}



#Replicates multistudysim_tot for 'reps' number of replications; all outputs represent aggregations of the outputs of multistudysim_tot across all reps
#Inputs: See earlier scripts
#coeffile: a string for the file name of the file containing the weights given to each tree  
#Outputs: 
#total_errors_tree, total_errors_forest: an array of the RMSE's or Log Losses for all weighting approaches across every iteration, for Weighting Trees and Weighting Forests respectively.
#total_coefs_tree, total_coefs_forest: all of the individual tree weights for all weighting approaches across every iteration, for Weighting Trees and Weighting Forests respectively. 
#means_tree, means_forest: column means of total_errors_tree and total_errors_forest, respectively 
#sds_tree, sds_forest: column sds of total_errors_tree and total_errors_forest, respectively 
#variance_inter_tree: array of the means and sds across all iterations of the percent variation in the outcome taken by interaction terms - this is equivalent for each iteration of  Weighting Trees and Weighting forests. 
#total_varimp_merged, total_varimp_tree, total_varimp_forest: the variable importance/tree structure tables for the Merged, Weighting Trees, and Weighting Forests approaches respectively. 
#For each level, total_varimp_tree and total_varimp_forest list the means over the top decile of tree weights across all variables, the means over the bottom decile, and the overall mean; these correspond to the three rows of the array for each level. 
rep_multistudy <- function(reps, modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, coeffile){
  ptm <- proc.time()
  total_errors_tree <- total_errors_forest <- matrix(NA, reps, 14)
  colnames(total_errors_tree) <- colnames(total_errors_forest) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                                                                    "Stack_noint", "Stack_noint_norm", "Stack_int",
                                                                    "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  
  total_varimp_tree <-  total_varimp_forest <-  matrix(NA, reps*3, 8)
  total_varimp_merged <-  matrix(NA, reps, 8)
  colnames(total_varimp_tree) <- colnames(total_varimp_forest) <- colnames(total_varimp_merged) <- c("tree_coefs", "num_var", "freq_truth", "num_inter", "prop_inter", "sum_varimp", "sum_totalimp", "sum_coef")
  for (i in 1:reps){
    ms_h <- multistudysim_helper(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc)
    mssi_tree <- c(multistudysim_tot(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight = TRUE, ms_h = ms_h))
    mssi_forest <- c(multistudysim_tot(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight = FALSE, ms_h = ms_h))
    
    total_errors_tree[i, ] <- mssi_tree$means
    total_errors_forest[i, ] <- mssi_forest$means
    
    j =  i + 2*(i-1)
    total_varimp_merged[i, ] <- ms_h$varimpsumm_mod0$means_tot
    total_varimp_tree[j:(j+2), ] <- rbind(mssi_tree$varimpsumm_mods$means_top10, mssi_tree$varimpsumm_mods$means_bottom10, mssi_tree$varimpsumm_mods$means_tot)
    total_varimp_forest[j:(j+2), ] <- rbind(mssi_forest$varimpsumm_mods$means_top10, mssi_forest$varimpsumm_mods$means_bottom10, mssi_forest$varimpsumm_mods$means_tot)
    if (i == 1){
      coefs_out_tree <- mssi_tree$coefs_list
      coefs_out_forest <- mssi_tree$coefs_list
      
      variance_inter_tree <- mssi_tree$variance_inter
      
    }
    else{
      coefs_out_tree <- rbind(coefs_out_tree, mssi_tree$coefs_list)
      coefs_out_forest <- rbind(coefs_out_forest, mssi_forest$coefs_list)
      
      variance_inter_tree <- rbind(variance_inter_tree, mssi_tree$variance_inter)
    }
    
    if (i %% 10 == 0){
      cat(paste0("iteration: ", i, "\n"))
    }
  }
  coefs_out_tree <- as.data.frame(coefs_out_tree)
  coefs_out_forest <- as.data.frame(coefs_out_forest)
  variance_inter_tree <- as.data.frame(variance_inter_tree)
  #variance_inter_forest <- as.data.frame(variance_inter_forest)
  
  coef_colnames <- c("coefs_stack_int", "coefs_stack_lasso", "coefs_stack_ridge", "coefs_ss_int", "coefs_ss_lasso", "coefs_ss_ridge")
  colnames(coefs_out_tree) <- colnames(coefs_out_forest) <- coef_colnames
  write.table(coefs_out_tree, paste0(coeffile,"_tree.csv"), sep = ",", row.names = F, col.names = coef_colnames)
  write.table(coefs_out_forest, paste0(coeffile,"_forest.csv"), sep = ",", row.names = F, col.names = coef_colnames)
  
  total_errors_tree <- na.omit(total_errors_tree)
  total_errors_forest <- na.omit(total_errors_forest)
  
  
  means_tree <- colMeans(total_errors_tree)
  means_forest <- colMeans(total_errors_forest)
  sds_tree <- apply(total_errors_tree, 2, sd)
  sds_forest <- apply(total_errors_forest, 2, sd)
  print(proc.time() - ptm)
  return(list(total_errors_tree = total_errors_tree, total_errors_forest = total_errors_forest, 
              total_coefs_tree = coefs_out_tree, total_coefs_forest = coefs_out_forest, 
              means_tree = means_tree, sds_tree = sds_tree,
              means_forest = means_forest, sds_forest = sds_forest,
              variance_inter_tree = variance_inter_tree,
              total_varimp_tree = total_varimp_tree, total_varimp_forest = total_varimp_forest,
              total_varimp_merged = total_varimp_merged)) 
}



#The overall simulation function
#runs rep_multistudy across varying levels of either heterogeneity or interaction strength
#inputs:
#var_list: a vector of coefficients to vary levels of (heterogeneity or interacction strength)
#het: TRUE if level heterogeneity being varied, FALSE if interaction strength being varied
#
#outputs:
#total_sds_tree, total_sds_forest: array of the sds across all iterations for each level in var_list, for Weighting Trees and Weighting Forests, respectively
#total_means_tree, total_means_forest: array of the means across all iterations for each level in var_list, for Weighting Trees and Weighting Forests, respectively
#total_variance_inter_tree: array of the means and sds across all iterations of the percent variation in the outcome taken by interaction terms, for every level - this is equivalent for each iteration of  Weighting Trees and Weighting forests. 
#For each level, the mean is the first row, the sd the second row; the columns represent the number of datasets with interactions 
#total_varimp_merged, total_varimp_tree, total_varimp_forest: the variable importance/tree structure tables for the Merged, Weighting Trees, and Weighting Forests approaches respectively. 
#For each level, and for every rep within a level, total_varimp_tree and total_varimp_forest list the means over the top decile of tree weights across all variables, the means over the bottom decile, and the overall mean; these correspond to the three rows of the array for each level. 
#For each level, total_varimp_merged lists the overall mean across all variables as well as the standard deviation. 
vary_levels <- function(reps, het, var_list, modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, allmeanstr, coefstr, meanstr, sdstr, varimpstr, varintstr){
  #errors_composite_tree <- errors_composite_forest <- array(0, c(reps, 14, length(var_list)))
  ptm = proc.time()
  
  total_means_tree <- total_means_forest <- array(0, c(length(var_list), 14))
  total_sds_tree <- total_sds_forest <- array(0, c(length(var_list), 14))
  total_variance_inter_tree <- array(0, c(2*length(var_list), ninter))
  total_varimp_merged <- array(0, c(2*length(var_list), 8))
  total_varimp_tree <- total_varimp_forest <- array(0, c(reps*3*length(var_list), 8))
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    if (het == TRUE){
      bad <- level
      val <- bad*.5}
    else{
      icoefs <- c(level, -level/3)
    }
    level_rep <- rep_multistudy(reps, modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc,  coeffile = paste0(coefstr, level, ".csv"))
    #errors_composite_tree[,,i] <- level_rep$total_errors_tree
    #errors_composite_forest[,,i] <- level_rep$total_errors_forest
    
    total_sds_tree[i, ] <- level_rep$sds_tree
    total_sds_forest[i, ] <- level_rep$sds_forest
    
    total_means_tree[i,] <- level_rep$means_tree
    total_means_forest[i,] <- level_rep$means_forest
    
    j = 2*i - 1
    total_variance_inter_tree[j, ] <- colMeans(as.matrix(level_rep$variance_inter_tree))
    total_variance_inter_tree[j+1, ] <- apply(as.matrix(level_rep$variance_inter_tree), 2, sd)
    
    total_varimp_merged[j, ] <- colMeans(as.matrix(level_rep$total_varimp_merged))
    total_varimp_merged[j+1, ] <- apply(as.matrix(level_rep$total_varimp_merged), 2, sd)
    
    k = (i-1)*(reps*3) + 1
    total_varimp_tree[k:(k + (reps*3 - 1)), ] <- as.matrix(level_rep$total_varimp_tree)
    total_varimp_forest[k:(k + (reps*3 - 1)), ] <- as.matrix(level_rep$total_varimp_forest)
    
  }
  colnames_totalmeans <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                           "Stack_noint", "Stack_noint_norm", "Stack_int",
                           "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  colnames(total_means_tree) <- colnames(total_means_forest) <- colnames(total_sds_tree) <- colnames(total_sds_forest) <- colnames_totalmeans
  rownames(total_means_tree) <- rownames(total_means_forest) <- rownames(total_sds_tree) <- rownames(total_sds_forest) <- var_list
  
  write.table(total_means_tree, paste0(meanstr,"_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_totalmeans)
  write.table(total_means_forest, paste0(meanstr,"_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_totalmeans)
  
  write.table(total_sds_tree, paste0(sdstr,"_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_totalmeans)
  write.table(total_sds_forest, paste0(sdstr,"_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_totalmeans)
  
  colnames_total_varimp <- c("tree_coefs", "num_var", "freq_truth", "num_inter", "prop_inter", "sum_varimp", "sum_totalimp", "sum_coef")
  write.table(total_varimp_tree, paste0(varimpstr,"_tree.csv"), sep = ",", row.names = rep(var_list, each = reps*3), col.names = colnames_total_varimp)
  write.table(total_varimp_forest, paste0(varimpstr,"_forest.csv"), sep = ",", row.names = rep(var_list, each = reps*3), col.names = colnames_total_varimp)
  write.table(total_varimp_merged, file = paste0(varimpstr,"_merged.csv"), sep = ",", row.names = rep(var_list, each = 2), col.names = colnames_total_varimp)
  
  colnames_totalvariance <- sapply(1:ninter, function(i) paste0("inter",i))
  if (ninter != 0){
    write.table(total_variance_inter_tree, paste0(varintstr,"_tree.csv"), sep = ",", row.names = rep(var_list, each = 2), col.names = colnames_totalvariance)
  }
  print(proc.time() - ptm)
  return(list(total_sds_tree = total_sds_tree, total_sds_forest = total_sds_forest, 
              total_means_tree = total_means_tree, total_means_forest = total_means_forest,
              total_variance_inter_tree = total_variance_inter_tree, 
              total_varimp_merged = total_varimp_merged, 
              total_varimp_tree = total_varimp_tree, total_varimp_forest = total_varimp_forest))
}

#Run the simulation for the three main population scenarios considered 

#vector of heterogeneity 'bad' level coefficients
het_ls <- seq(0, 10, .5)

#vector of interaction strength coefficients
inter_ls <- seq(2, 12, .5)

#(a) 2 interactions in training, 2 in testing
case1_scen1 <- vary_levels(reps = 100, het = FALSE, var_list = ntree_ls, modfit = randomforestfit, modpred = randomforestpredict, good = .25, bad = 1, val = .4, edat_orig, simtype = "nonl", complexity = 10, numtree = 10, bin = TRUE, ninter = 2, icoefs = c(4.4, -1.8), setnc = TRUE, allmeanstr = "allmeans_c1s2_", coefstr = "coefs_c1s2_", meanstr = "means_c1s2",  zerosstr = "zeros_c1s2", sdstr = "sds_c1s2")

#(b) 6 interactions in training, 2 in testing
case1_scen2 <- vary_levels(reps = 100, het = FALSE, var_list = inter_ls, modfit = randomforestfit, modpred = randomforestpredict, good = .25, bad = 1, val = .4, edat_orig, simtype = "nonl", complexity = 10, numtree = 10, bin = TRUE, ninter = 8, icoefs = c(4.4, -1.8), setnc = TRUE, allmeanstr = "allmeans_c1s3_", coefstr = "coefs_c1s3_", meanstr = "means_c1s3",  zerosstr = "zeros_c1s3", sdstr = "sds_c1s3")

#(c) Varying level of heterogeneity, no interactions: 
case1_scen3 <- vary_levels(reps = 100, het = TRUE, var_list = het_ls, modfit = randomforestfit, modpred = randomforestpredict, good = .25, bad = 1, val = .4, edat_orig, simtype = "nonl", complexity = 10, numtree = 10, bin = TRUE, ninter = 0, icoefs = c(4.4, -1.8), setnc = TRUE, allmeanstr = "allmeans_c1s4_", coefstr = "coefs_c1s4_", meanstr = "means_c1s4",  zerosstr = "zeros_c1s4", sdstr = "sds_c1s4")


#These simulations are used to create Figure 2, Table 1, 
