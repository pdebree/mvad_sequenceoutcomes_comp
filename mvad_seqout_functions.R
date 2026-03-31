# Functions for Sequence Outcome Prediction Work - CFDA Competition, Linear and Non-Linear

create_dists <- function(data.seq) {
  # Creates a substition cost matrix, using Transition Rate (T-Rate) -  which gets the substitution costs from 
  # the transition probabilities from all pairs of states (higher probability means
  # lower cost). Here we just use it to get the size of a substitution matrix quickly 
  
  # Calculate costs 
  costs.data <- seqcost(data.seq, method = "INDELSLOG", with.missing = TRUE)

  submat   <- seqsubm(data.seq, method= "TRATE")
  submat.2 <- submat
  submat.2[] <- 1
  diag(submat.2) <- 0
  
  # OM - T-RATE - Simply use the transition rates between each state. 
  dist.om_trate <- seqdist(data.seq, method="OM", indel = 1, sm = submat, with.missing = TRUE)
  # LCS
  dist.lcs <- seqdist(data.seq, method="LCS")
  # OM - SLOG 
  dist.om_slog <- seqdist(data.seq, method="OM", indel = costs.data$indel, sm = costs.data$sm, with.missing = TRUE)

  dists <- list(`OM-trate`=dist.om_trate,`LCS`=dist.lcs, `OM-slog`=dist.om_slog)
  
  return(dists)
}




# Assigns a hard cluster based on a sequence (with similarity matrix) 
assign_new <- function(dist_mat,train_idx,train.clust,test_idx) {
  
  #spread out the clusters and assign fake label to test set
  len <- length(train_idx) #boolean vector
  clusts <- rep(0,len) 
  
  clusts[train_idx] <- train.clust # assigned labels (with 0 being the new data)
  
  dist2clust <- apply(dist_mat,1,function(x,cl) {tapply(x,cl,mean)},cl=clusts)  
  
  newClusts <- apply(dist2clust[-1,],2,which.min) #drop first row (unassigned cluster)
  
  return(newClusts[test_idx])
  
}

# Assigns soft cluster membership probabilities 
assign_new_fanny <- function(dist_mat,train.idx,train.memb,test.idx,memb.exp=1.5,n.draws=300,prior.n.equiv=100) {
  
  len <- length(train.idx) # all members 
  simmat <- 1-dist_mat/max(dist_mat) # all similarities 
  test_locs <- (1:len)[test.idx] # indices of test locations
  
  train.memb <- as.matrix(train.memb) # train membership (in cluster) probabilities 
  K <- ncol(train.memb) # number of clusters 
  M <- length(test_locs) # number of new obs to assign
  test.memb <- matrix(NA,M,K) # empty membership probability matrix 
  
  for (i in seq_along(test_locs)) {
    
    sim.i <- simmat[test_locs[i],train.idx] # similarities of ith test row with all training data points
    sim.i <- sim.i/sum(sim.i) # relative similarity (along row)
    wtd.memb <- sim.i%*%train.memb # essentially "weighing" test similarities by training cluster memberships 
    
    # pull dirichlet draws based on weighted memberships to get probabilites 
    u.prop <- gtools::rdirichlet(n.draws,prior.n.equiv*wtd.memb) # n.draws - make 200
    d.i <- dist_mat[test_locs[i],train.idx,drop=F] # distance values to each train point
    
    numer <- d.i%*%(train.memb^memb.exp)%*%t(u.prop^memb.exp) 
    denom <- apply(u.prop^memb.exp,1,sum) #only relevant terms
    
    # use draw that minimizes 
    min.loc <- which.min(numer/denom)
    test.memb[i,] <- u.prop[min.loc,]
    
  }
  return(test.memb)
}


# soft clust assign probs based on medoids (FKM.med):
pred.fclust <- function(dist2med,m=1.5) {
  
  # Handle the 'Distance = 0' edge case 
  eps <- 1e-10
  dist2med[dist2med==0] <- eps
  
  # Calculate the weights (1 / distance ^ power)
  # The power is 2/(m-1)
  power <- 2 / (m - 1)
  weights <- (1 / dist2med)^power
  
  # Normalize row-wise so probabilities sum to 1
  # This gives you the soft membership matrix
  test_probs <- weights / rowSums(weights)
  return(test_probs)
}




train_test_lm_clust_hard <- function(clustering, dists, nClusts, dat, Y.cont, train_idx,test_idx) {
  cut1 <- cutree(clustering,k=nClusts)
  clust1.fac <- factor(cut1)
  train_dat.cont <- dat[train_idx,] %>% mutate(cluster=clust1.fac, Y=Y.cont[train_idx])
  
  if (nClusts==1) {
    cut2 <- rep(1,sum(test_idx)) #the trivial cluster
  } else cut2 <- assign_new(dists,train_idx,cut1,test_idx)
  
  clust2.fac <- factor(cut2)
  test_dat.cont <- dat[test_idx,] %>% mutate(cluster=clust2.fac, Y=Y.cont[test_idx])
  
  if (nClusts==1) {
    # remove clusters 
    fit.lm <- lm(Y~., data=train_dat.cont%>%dplyr::select(-cluster)) 
  } else {
    fit.lm <- lm(Y~., data=train_dat.cont) 
  }
  
  preds.cont <- predict(fit.lm, newdata = test_dat.cont) 
  mse <- sum((preds.cont - Y.cont[test_idx])^2)/nrow(test_dat.cont)
  
  return(mse)
  
}

train_test_lm_clust_soft <- function(dist_matrix, nClusts, dat, Y.cont, train_idx,test_idx, fuzziness=1.5) {
  
  
  clustering <- fanny(dist_matrix[train_idx, train_idx], k=nClusts, memb.exp=fuzziness, diss=TRUE,maxit = 1000)
  
  train.memb <- clustering$membership
  test.memb <- assign_new_fanny(dist_mat=dist_matrix,train.idx=train_idx, train.memb=train.memb, 
                               test.idx=test_idx, memb.exp=fuzziness,n.draws=200,prior.n.equiv=100) 
  
  if (nClusts == 2) {
    train_dat.cont <- dat[train_idx,] %>% mutate(Y=Y.cont[train_idx], cluster1=train.memb[,-1])
    test_dat.cont <- dat[test_idx,] %>% mutate(Y=Y.cont[test_idx], cluster1=test.memb[,-1])
  } else {
    train_dat.cont <- cbind(dat[train_idx,], train.memb[,-1]) %>% mutate(Y=Y.cont[train_idx])
    test_dat.cont <- cbind(dat[test_idx,], test.memb[,-1]) %>% mutate(Y=Y.cont[test_idx])
  } 
  
  fit.lm <- lm(Y~., data=train_dat.cont)
  
  preds.cont <- predict(fit.lm, newdata = test_dat.cont)
  mse <- sum((preds.cont - Y.cont[test_idx])^2)/nrow(test_dat.cont)
  
  
  preds.train <- predict(fit.lm)
  
  output <- list()
  output$mse <- mse
  output$converged <- clustering$convergence["converged"] == 1

  return(output)
  
}

train_test_lm_clust_soft_fkmmed <- function(dist_matrix, nClusts, dat, Y.cont, train_idx,test_idx, m=1.5, RS=20) {
  # default fits has m=1.5m RS=20 (k is number of clusters - passed through the nClusts parameter.)
  
  # fit FKM.med with the training data 
   
  fcl <- FKM.med(as.dist(dist_matrix[train_idx, train_idx]),k=nClusts,m=m,RS=RS)
  

  train.memb <- fcl$U 
  
  # tried this - different outcomes - not sure why  
  train.memb <- pred.fclust(dist_matrix[train_idx,fcl$medoid])
  
  # pull medoids from fit to "predict" test set cluster memberships (based on training medoids)
  test.memb <- pred.fclust(dist_matrix[test_idx,fcl$medoid])
  
  # update names of test.memb columns to match train.memb (for prediction modeling internals)
  colnames(test.memb) <- colnames(train.memb)
  
  if (nClusts == 2) {
    train_dat.cont <- dat[train_idx,] %>% mutate(Y=Y.cont[train_idx], cluster1=train.memb[,-1])
    test_dat.cont <- dat[test_idx,] %>% mutate(Y=Y.cont[test_idx], cluster1=test.memb[,-1])
  } else {
    train_dat.cont <- cbind(dat[train_idx,], train.memb[,-1]) %>% mutate(Y=Y.cont[train_idx])
    test_dat.cont <- cbind(dat[test_idx,], test.memb[,-1]) %>% mutate(Y=Y.cont[test_idx])
  } 
  
  fit.lm <- lm(Y~., data=train_dat.cont)
  # train preds
  preds.train <- predict(fit.lm)
  
  # predict test data
  preds.cont <- predict(fit.lm, newdata = test_dat.cont)
  
  mse <- sum((preds.cont - Y.cont[test_idx])^2)/nrow(test_dat.cont)
  
  return(mse)

}




which.min.arr <- function(x) arrayInd(which.min(x), dim(x))

train_test_lm_comps <- function(dat, Y.cont, Y.bin, train_idx, test_idx, end_col) {
  
  # This is fairly hard coded to the form of mvad_covars0. 
  train_dat.cont <- dat[train_idx,c(1:end_col)] %>% mutate(Y=Y.cont[train_idx])
  test_dat.cont <- dat[test_idx,c(1:end_col)] %>% mutate(Y=Y.cont[test_idx])
  
  fit.lm <- lm(Y~., data=train_dat.cont)
  
  preds.cont <- predict(fit.lm, newdata = test_dat.cont) 
  mse <- sum((preds.cont - Y.cont[test_idx])^2)/nrow(test_dat.cont)
  
  list(mse=mse,fit.lm=fit.lm)
}


train_test_lm_comps_mets <- function(train_dat, test_dat, Y.cont, Y.bin, train_idx, test_idx, end_col) {
  
  # This is fairly hard coded to the form of mvad_covars0. 
  train_dat.cont <- train_dat[,1:end_col] %>% mutate(Y=Y.cont[train_idx])
  test_dat.cont <- test_dat[,1:end_col] %>% mutate(Y=Y.cont[test_idx])
  
  fit.lm <- lm(Y~., data=train_dat.cont)
  
  preds.cont <- predict(fit.lm, newdata = test_dat.cont) 
  mse <- sum((preds.cont - Y.cont[test_idx])^2)/nrow(test_dat.cont)
  
  list(mse=mse,fit.lm=fit.lm)
}

train_test_lm_clust_hard_seq <- function(clustering, dist_matrix, nClusts, dat, Y.cont, train_idx,test_idx,train_seq_pcs, test_seq_pcs) {
  cut1 <- cutree(clustering,k=nClusts)
  clust1.fac <- factor(cut1)
  
  train_covars <- cbind(dat[train_idx,], train_seq_pcs)
  train_dat.cont <- train_covars %>% mutate(cluster=clust1.fac, Y=Y.cont[train_idx])
  
  if (nClusts==1) {
    cut2 <- rep(1,sum(test_idx)) #the trivial cluster
  } else cut2 <- assign_new(dist_matrix,train_idx,cut1,test_idx)
  
  clust2.fac <- factor(cut2)
  
  # add in clusters - pcs as covars 
  test_covars <- cbind(dat[test_idx,], test_seq_pcs)
  test_dat.cont <- test_covars %>% mutate(cluster=clust2.fac, Y=Y.cont[test_idx])
  
  if (nClusts==1) {
    fit.lm <- lm(Y~., data=train_dat.cont%>%dplyr::select(-cluster)) 
  } else {
    fit.lm <- lm(Y~., data=train_dat.cont) 
  }
  
  preds.cont <- predict(fit.lm, newdata = test_dat.cont) 
  mse <- sum((preds.cont - Y.cont[test_idx])^2)/nrow(test_dat.cont)
  list(mse=mse,fit.lm=fit.lm)
}

train_test_lm_clust_soft_seq <- function(dist_matrix, nClusts, dat, Y.cont, train_idx,test_idx, fuzziness=2, train_seq_pcs, test_seq_pcs) {
  
  clustering <- fanny(dist_matrix[train_idx, train_idx], k=nClusts, memb.exp=fuzziness, diss=TRUE,maxit = 1000)
  
  train.memb <- clustering$membership
  
  test.memb <- assign_new_fanny(dist_mat=dist_matrix,train.idx=train_idx, train.memb=train.memb,
                               test.idx=test_idx, memb.exp=fuzziness,n.draws=200,prior.n.equiv=100) 
  
  if (nClusts == 2) {
    train_dat.cont <- cbind(dat[train_idx,], train_seq_pcs) %>% mutate(Y=Y.cont[train_idx], cluster1=train.memb[,-1])
    test_dat.cont <- cbind(dat[test_idx,], test_seq_pcs) %>% mutate(Y=Y.cont[test_idx], cluster1=test.memb[,-1])
  } else {
    train_dat.cont <- cbind(dat[train_idx,], train_seq_pcs, train.memb[,-1]) %>% mutate(Y=Y.cont[train_idx])
    test_dat.cont <- cbind(dat[test_idx,], test_seq_pcs, test.memb[,-1]) %>% mutate(Y=Y.cont[test_idx])
  } 
  
  fit.lm <- lm(Y~., data=train_dat.cont)
  
  preds.cont <- predict(fit.lm, newdata = test_dat.cont)
  mse <- sum((preds.cont - Y.cont[test_idx])^2)/nrow(test_dat.cont)
  
  return(list(mse=mse,fit.lm=fit.lm))
}


train_mse_rf <- function(train_data, test_data, mtry = NA) {

  if (!is.na(mtry)) {
    mtry_tr = mtry
  } else {
    mtry_tr = ncol(train_data) - 1
  }
 
  # fit with hard clusters 
  fit.rf <- ranger(num_month_em_last_year ~ .,
                        data = train_data,
                        num.trees = 1000, 
                        mtry = mtry_tr,
                        respect.unordered.factors = TRUE)


  rf_pred <- predict(fit.rf, data = test_data)
  
  return(sum((preds <- rf_pred$predictions - test_data$num_month_em_last_year)^2)/ nrow(test_data))
  
}


hard_cluster <- function(clusterward, nClusts, covars, num_month_em_last_year, train_idx, test_idx, dist_matrix) {
  
  cut1 <- cutree(clusterward,k=nClusts)
  clust1.fac <- factor(cut1)
  train_data <- covars[train_idx,] %>% mutate(cluster=clust1.fac, num_month_em_last_year=num_month_em_last_year[train_idx])
  
  cut2 <- assign_new(dist_matrix,train_idx,cut1,test_idx)
  clust2.fac <- factor(cut2)
  test_data <- covars[test_idx,] %>% mutate(cluster=clust2.fac, num_month_em_last_year=num_month_em_last_year[test_idx])
  
  return(list(train_data=train_data, test_data=test_data))
}
  
soft_cluster <- function(dist_matrix,train_idx,test_idx, nClusts, fuzziness=1.5, covars, num_month_em_last_year) {
  
  clustering_soft <- fanny(dist_matrix[train_idx, train_idx], 
                           k=nClusts, memb.exp=fuzziness, diss=TRUE, maxit = 1000)
  
  train.memb <- clustering_soft$membership
  test.memb <- assign_new_fanny(dist_mat=dist_matrix,train.idx=train_idx, train.memb=train.memb, 
                               test.idx=test_idx,memb.exp=fuzziness,n.draws=200,
                               prior.n.equiv=100) 
  
  colnames(train.memb) <- paste0("Cluster",1:nClusts)
  colnames(test.memb) <- paste0("Cluster",1:nClusts)
  
  if (nClusts == 2) {
    train_data <- covars[train_idx,] %>% 
      mutate(num_month_em_last_year=num_month_em_last_year[train_idx], cluster1=train.memb[,-1])
    test_data <- covars[test_idx,] %>% 
      mutate(num_month_em_last_year=num_month_em_last_year[test_idx], cluster1=test.memb[,-1])
  } else {
    train_data <- cbind(covars[train_idx,], train.memb[,-1]) %>%
      mutate(num_month_em_last_year=num_month_em_last_year[train_idx])
    test_data <- cbind(covars[test_idx,], test.memb[,-1]) %>%
      mutate(num_month_em_last_year=num_month_em_last_year[test_idx])
  } 
  
  return(list(train_data=train_data, test_data=test_data))
}


pivot_means <- function(fold_means, comp_name, nComps) {
  dimnames(fold_means) <- list(1:nComps, 1:11)
  
  fold_means_long <- fold_means %>%
    as.data.frame() %>%
    mutate(index = row_number()) %>%             
    pivot_longer(cols = -index, 
                 names_to = "Mtry", 
                 values_to = "RMSE") %>% rename({{comp_name}}:=index) %>% arrange(Mtry)
  
  fold_means_long$Mtry <- as.factor(fold_means_long$Mtry)
  return(fold_means_long)
}

plot_comps_mtry <- function(fold_means, comp_lab){
  ggplot(data=fold_means, aes(x=.data[[names(fold_means)[1]]], y=RMSE, color=as.factor(Mtry))) + 
    geom_line() + 
    labs(title=paste0("Performance with Varying Number of ", comp_lab), 
         color="Mtry", x=comp_lab) + 
    scale_color_discrete(breaks = seq(1,11, 1))
  
}

calc_mse <- function(preds, y_dat) {
  return(sum((preds - y_dat)^2)/length(y_dat))
  
}

# CHANGE THIS
fit_linear_train_test <- function(train_data, test_data) {
  fit.lin <- lm(y~., data=train_data) 
  preds.lin <- predict(fit.lin, newdata = test_data) 
  mse.lin_test <- calc_mse(preds.lin, y_dat=test_data[,"y"]) 
  
  preds.lin_train <- predict(fit.lin, newdata = train_data) 
  mse.lin_train <- calc_mse(preds.lin_train, y_dat=train_data[,"y"]) 
    
  
  output <- list()
  output$train <- mse.lin_train
  output$test <- mse.lin_test
  return(output)
  
}

fit_linear <- function(train_data, test_data) {
  fit.lin <- lm(y~., data=train_data) 
  preds.lin <- predict(fit.lin, newdata = test_data) 
  mse.lin <- calc_mse(preds.lin, y_dat=test_data[,"y"]) 

  return(mse.lin)
  
}







