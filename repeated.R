#create datasets by repeating the covariate values for one 10 times
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

#Requires to have first saved the ovarian cancer datasets: in this case, it is saved as edat_orig.Rda
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

init_data <- function(edat_orig, ndat, nvar, simtype, ninter, good, bad, val, icoefs, setnc, bin){
  
  edat <- edat_orig
  edat <- edat[sample(1:ndat)] 
  
  idx <- sample(1:ncol(edat[[1]]), nvar)
  repeated_dataset <- edat[[1]][,idx]
  for(i in 1:10){
    edat[[i]] <- repeated_dataset
  }
  for(i in 11:ndat){
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
  
  
  for(i in 1:ndat){ #ndat is between 1 and 15 (number of datasets)
    
    if (bad == 0){
      curcoefs <- coefs
    } else{
      curcoefs <- sapply(coefs, function(x){runif(1, x - (i<=5)*good - 
                                                    (i > 5 & i <= 10)*bad - (i > 10)*val , x + 
                                                    (i<=5)*good + (i > 5 & i <= 10)*bad + (i > 10)*val)}) #adds noise to the coefficients
    }
    
    if(simtype == "slash"){
      y <- (edat[[i]][,vars] %*% curcoefs) + 
        0.05*cbind(rnorm(nrow(edat[[i]]))/runif(nrow(edat[[i]]))) # Added "slash" noise
    } else if(simtype == "nonl"){
      if ((1 <= i) & (i <= ninter)){
        y <- (edat[[i]][,vars] %*% curcoefs) + icoefs[1]*edat[[i]][,vars[1]]*edat[[i]][,vars[2]] 
        - icoefs[2]*edat[[i]][,vars[1]]*edat[[i]][,vars[3]] + cbind(rnorm(nrow(edat[[i]]))) # Added interaction terms
      }
      else{
        y <- (edat[[i]][,vars] %*% curcoefs) + cbind(rnorm(nrow(edat[[i]]))) # Added noise
      }
    } else {
      y <- (edat[[i]][,vars] %*% curcoefs)
      #y <- (edat[[i]][,vars] %*% curcoefs) + cbind(rnorm(nrow(edat[[i]]))) # Added noise
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
  edat
}

#complexity = max depth allowed in each decision tree
#numtree = number of trees in each forest
#bin = TRUE if using binary outcome data, false if it is continuous
#ninter: number of datasets out of the 15 with interactinos present: floor(ninter/3) is the number of testing datasets with interactions, and ninter - floor(ninter/3) is the number of training datasets with interactions. 
#icoefs: a vector of two coefficients corresponding to two interactions between variables present in init_data if ninter > 1
#setnc = TRUE if setting the number of coefficients used to generate the outcome to a particular number in init_data, false if it is a randomly chosen number between 2/3 and nvar
#treeweight = TRUE if weighting individual trees, false if weighting forests

multistudysim_tot <- function(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight){
  
  ndat <- length(edat_orig)
  ntrain <- 10
  nvar <- 100
  modfit <- modfit
  modpred <- modpred
  
  edat <- init_data(edat_orig, ndat, nvar, simtype, ninter, good, bad, val, icoefs, setnc, bin)
  
  if (ninter == 2){ #one study with interaction in training, one in testing
    edat <- replace(edat, c(2, 11), edat[c(11, 2)])
  } else if (ninter > 2){
    numrepl <- floor(ninter/3)
    edat <- replace(edat, (ninter - numrepl):ninter, edat[11:(11 + numrepl)])
  }
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
    merged <- modpred(mod0, newdata = edat[[i]][,-1], treeweight = FALSE)	
    merged <- as.vector(sapply(merged, as.numeric))
    if (treeweight == TRUE){
      #merged <- apply(merged, 2, as.numeric)
      allmod <- t(do.call(cbind,lapply(mods, modpred, newdata = edat[[i]][,-1], treeweight)))
      allmod <- apply(allmod, 2, as.numeric)
    }
    else{
      #merged <- as.vector(sapply(merged, as.numeric))
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

rep_multistudy <- function(reps, modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight, allmeansfile, coeffile){
  ptm <- proc.time()
  total_errors <- matrix(NA, reps, 13)
  zeros_summary <- matrix(NA, reps, 6)
  colnames_totalerror <- c("Merged", "Unweighted", "CS_Weighted",
                           "Stack_noint", "Stack_noint_norm", "Stack_int",
                           "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  colnames(total_errors) <- colnames_totalerror
  colnames(zeros_summary) <- c("stack_0", "stacklasso_0", "SS_0", "SSlasso_0", "stack_intersect", "SS_intersect")
  
  for (i in 1:reps){
    mssi <- c(multistudysim_tot(modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight))
    total_errors[i, ] <- mssi$means[-3]
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
  print(proc.time() - ptm)
  return(list(total_errors = total_errors, zeros = zeros_summary, total_coefs = coefs_out, means = means, sds = sds)) 
}


#var_list: a vector of coefficients to vary levels of (heterogeneity or interacction strength)
#het = TRUE if level heterogeneity being varied, FALSE if interaction strength being varied
vary_levels <- function(reps, het, var_list, modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight, allmeanstr, coefstr, meanstr, zerosstr, sdstr){
  errors_composite <- array(0, c(reps, 13, length(var_list)))
  total_zeros <- array(0, c(length(var_list),  6))
  total_means <- array(0, c(length(var_list), 13))
  total_sds <- array(0, c(length(var_list), 13))
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    if (het == TRUE){
      bad <- level
      val <- bad*.5}
    else{
      icoefs <- c(level, -level/3)
    }
    level_rep <- rep_multistudy(reps, modfit, modpred, good, bad, val, edat_orig, simtype, complexity, numtree, bin, ninter, icoefs, setnc, treeweight, allmeansfile = paste0(allmeanstr, level, ".csv"), coeffile = paste0(coefstr, level, ".csv"))
    
    errors_composite[,,i] <- level_rep$total_errors
    total_sds[i, ] <- level_rep$sds
    total_zeros[i, ] <- level_rep$zeros
    total_means[i,] <- level_rep$means
    print("means")
    print(level_rep$means)
    print("sds")
    print(level_rep$sds)
    print("zeros")
    print(level_rep$zeros)
  }
  colnames_totalmeans <- c("Merged", "Unweighted", "CS_Weighted",
                           "Stack_noint", "Stack_noint_norm", "Stack_int",
                           "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  colnames_zeros <- c("stack_0", "stacklasso_0", "SS_0", "SSlasso_0", "stack_intersect", "SS_intersect")
  colnames(total_means) <- colnames_totalmeans
  colnames(total_sds) <- colnames_totalmeans
  colnames(total_zeros) <- colnames_zeros
  rownames(total_means) <- var_list
  rownames(total_sds) <- var_list
  rownames(total_zeros) <- var_list
  write.table(total_means, meanstr, sep = ",", row.names = var_list, col.names = colnames_totalmeans)
  write.table(total_sds, sdstr, sep = ",", row.names = var_list, col.names = colnames_totalmeans)
  write.table(total_zeros, zerosstr, sep = ",", row.names = var_list, col.names = colnames_zeros)
  return(list(errors_composite = errors_composite, total_sds = total_sds, total_zeros = total_zeros, total_means = total_means))
}


#vector of heterogeneity 'bad' level coefficients
het_ls <- seq(0, 10, .5)

##Varying level of heterogeneity: 
### Scenario 4: no interactions
case1_scen4 <- vary_levels(reps = 100, het = TRUE, var_list = het_ls, modfit = randomforestfit, modpred = randomforestpredict, good = 0, bad = 0, val = 0, edat_orig, simtype = "normal", complexity = 10, numtree = 10, bin = FALSE, ninter = 0, icoefs = c(4.4, -1.8), setnc = TRUE, treeweight = TRUE, allmeanstr = "allmeans_c1s4_", coefstr = "coefs_c1s4_", meanstr = "means_c1s4",  zerosstr = "zeros_c1s4", sdstr = "sds_c1s4")

