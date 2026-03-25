# CFDA Linear Competition 
# First written for Multi-Level Models practicum in Fall 2025

library(tidyverse)
library(TraMineR)
library(cluster)
library(cfda)
source("mvad_seqout_functions.R")


### Data Preparation 
# Data load in and labeling 
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
                 "joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")

# seqdef - creates a state sequence object (formulates sequences into labels)
mvad.seq <- seqdef(mvad, 15:50, states=mvad.scodes, labels=mvad.labels)

# Calculate number of months employed in last year. 
mvad_last_year <- mvad[,75:86]
num_month_em_last_year <- apply(mvad_last_year, 1, function(x) length(which(x=="employment")))
mvad_covars <- mvad[3:14] %>% dplyr::select(-Western) #reference group

# Create distance matrices 
dists <- create_dists(data.seq=mvad.seq)



### Training set up 
folds <-  5 
nrecs <- nrow(mvad[1])

# Set up training and test split so it is the same across the competition 
idx.orig <- rep(1:folds,each=floor(nrecs/folds))
if (nrecs %% folds != 0) idx.orig <- c(idx.orig,1:(nrecs %% folds))
idx <- sample(idx.orig) #shuffle

# Creates a random shuffle of the indices to be used in the the train-test split and folds. 
shuffle_index <- sample(1:nrecs, nrecs)
ids <- sort(unique(mvad$id))

# Number of each method to check
nClusts <- 25 
nWindows <- 16 
nHarms <- 25

### Arrays to hold RMSEs from different runs 
rmse.cv.om_hard <- array(NA,c(folds,nClusts))
rmse.cv.om_soft <- array(NA,c(folds,nClusts))
rmse.cv.lcs_hard <- array(NA,c(folds,nClusts))
rmse.cv.lcs_soft <- array(NA,c(folds,nClusts))
rmse.cv.windows <- array(NA,c(folds,nWindows))
rmse.cv.om_seq_hard <- array(NA,c(folds,nClusts))
rmse.cv.om_seq_soft <- array(NA,c(folds,nClusts))
rmse.cv.harm <- array(NA,c(folds,nHarms))




### Cross Validation for Clustering Methods with OM-Transition Rate 
for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx
  
  # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
  clusterward_hard <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")
  
  # hard coded for OM-Trate
  for (j in 2:nClusts) {
    
    # hard clustering
    rmse.cv.om_hard[i,j] <- train_test_lm_clust_hard(clustering=clusterward_hard, dists=dists[[1]], nClusts=j,  dat=mvad_covars, Y.cont=num_month_em_last_year, train_idx=train_idx, test_idx=test_idx)
    
    #soft clustering - fanny 
    if (j < 12) {
      rmse.cv.om_soft[i,j] <- train_test_lm_clust_soft(dist_matrix=dists[[1]], nClusts=j, dat=mvad_covars, Y.cont=num_month_em_last_year, train_idx=train_idx, test_idx=test_idx, fuzziness = 1.5)
    }
  }
}


### Cross Validation for Clustering Methods with LCS as a Distance Measure
for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx
  
  # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR LCS-DSS
  clusterward <- agnes(dists[[2]][train_idx,train_idx], diss=TRUE, method="ward")
  
  for (j in 2:nClusts) {
    rmse.cv.lcs_hard[i,j] <- train_test_lm_clust_hard(clustering=clusterward, dists=dists[[2]], 
                                                      nClusts=j,  dat=mvad_covars, Y.cont=num_month_em_last_year,
                                                      train_idx=train_idx, test_idx=test_idx)
    
    # soft clustering - removed because of convergence errors
    if (j < 15) {
      rmse.cv.lcs_soft[i,j] <- train_test_lm_clust_soft(dist_matrix = dists[[2]], nClusts=j, dat=mvad_covars,
                                                        Y.cont=num_month_em_last_year,
                                                        train_idx=train_idx,
                                                        test_idx=test_idx, fuzziness = 1.25)
    }
  }
}

### Year-State Count Windows (PCA)

# Make our windows of data - first three years
mvad_states <- mvad[,c(1,15:50)]

# column names for the three years 
y1 <- colnames(mvad_states)[2:13]
y2 <- colnames(mvad_states)[14:25]
y3 <- colnames(mvad_states)[26:37]

# makes year and id combinations with counts for each state
id_year <- mvad_states %>% pivot_longer(cols=colnames(mvad_states[2:37]), names_to="output") %>% 
  mutate(year = case_when(output %in% y1 ~ 1, output %in% y2 ~ 2, TRUE ~ 3)) %>% count(id, year, value)

# make full table with 16 columns (one for each year/id combo) + 1 for id.
mvad_states_wide <- id_year %>% 
  pivot_wider(id_cols=c(id), 
              names_from=c(value,year), values_from=n, values_fill = 0) %>% dplyr::select(-id)

mvad_windows <- cbind(mvad_covars, mvad_states_wide)


for (i in 1:folds) {

  cat("Fold Number ",i,"\n")
  test_idx <- idx == i  
  train_idx <- !test_idx
  
  # pull out training and testing for only the windows
  train_windows <- mvad_windows[12:27][train_idx, ]
  test_windows <- mvad_windows[12:27][test_idx, ]
  
  pca_windows_train <- prcomp(x=train_windows, center=TRUE, scale=TRUE)
  train_scores <- pca_windows_train$x
  test_scores <- predict(pca_windows_train, newdata = test_windows) 
  
  # add principal components to demographic data (no longer need original windows)
  train_mvad_windows <- cbind(mvad_windows[1:11][train_idx,], train_scores)
  test_mvad_windows <- cbind(mvad_windows[1:11][test_idx,], test_scores)
  
  for (j in 1:nWindows) {
    # hard coded for the format of the mvad data 
    # j - 1 because we do want to look at only the 
    # based on the variance decreasing through the pcas, we can go through the 
    # components based on index. 
    # appropriate training and testing data is pull out based on the 
    rslt <- train_test_lm_comps(
      dat=mvad_windows, Y.cont=num_month_em_last_year, 
      train_idx=train_idx, test_idx=test_idx, end_col = 11+j) 
    rmse.cv.windows[i,j] <- rslt$rmse
  }
}


### Sequence PCs + Clustering Methods 

mvad_states <- mvad[,15:50]

# encodings for the statebadness 
st_alphabet <- alphabet(mvad.seq) 
# Example: If alphabet is "A", "B", "C"
st_prec_values <- c(1, -1, -1, 2, -1, -1) # A=1 (low badness), C=3 (high badness)
names(st_prec_values) <- st_alphabet

mvad_rmetrics <- mvad_states %>% mutate(
  spells=seqindic(mvad.seq, "dlgth")$Dlgth, 
  visited_states=seqindic(mvad.seq, "visited")$Visited, 
  num_of_trans = seqindic(mvad.seq,"trans")$Trans, 
  mean_spell_dur = seqindic(mvad.seq,"meand")$MeanD, 
  # pedantic - this one pulled out a "seqivardur" "numeric" datatype (not sure 
  # how it did both, so I have to force it to be numeric)
  sd_spell_dur = as.numeric(seqindic(mvad.seq,"dustd")$Dustd), 
  # Diversity I
  entropy = seqindic(mvad.seq,"entr")$Entr, 
  # more interested in states than spells (see paper)
  dss_subs = seqindic(mvad.seq,"nsubs")$Nsubs, 
  complexity = seqindic(mvad.seq,"cplx")$Cplx, 
  # could look at other turbulence measures 
  turbulence = seqindic(mvad.seq,"turb")$Turb, 
  badness = seqibad(seqdata = mvad.seq,stprec = st_prec_values), 
  degradation = seqidegrad(seqdata = mvad.seq, stprec = st_prec_values), 
  insecurity = seqinsecurity(seqdata = mvad.seq, stprec = st_prec_values))


rmetrics <- mvad_rmetrics[,37:48] 
nSeqPcs <- ncol(rmetrics) - 1


# go through folds in repeat
for (i in 1:folds) {
  cat("Fold Number ",i,"\n")
  test_idx <- idx == i
  train_idx <- !test_idx
  
  # using rmetrics from above
  train_rmetrics <- rmetrics[train_idx, ]
  test_rmetrics <- rmetrics[test_idx, ]
  
  pca_comps_train <- prcomp(x=train_rmetrics, center=TRUE, scale=TRUE)
  
  train_scores <- pca_comps_train$x
  test_scores <- predict(pca_comps_train, newdata = test_rmetrics) 
  
  
  # create a agnes tree for the clusters, based on the 1st similarity matrix - HARD CODED FOR OM
  clusterward_hard <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")
  
  
  # keep as matrices, for easier combination with covariate data 
  train_k_scores <- train_scores[, 1:nSeqPcs, drop = FALSE]
  test_k_scores <- test_scores[, 1:nSeqPcs, drop = FALSE]
  
  for (j in 2:nClusts) {
    
    # hard clustering
    rslt_hard <- train_test_lm_clust_hard_seq(clustering=clusterward_hard, dist_matrix=dists[[1]], nClusts=j,  dat=mvad_covars, Y.cont=num_month_em_last_year,
                                              train_idx=train_idx, test_idx=test_idx, 
                                              train_seq_pcs=train_k_scores, test_seq_pcs = test_k_scores)
    rmse.cv.om_seq_hard[i,j] <- rslt_hard$rmse
    
    #soft clustering
    if (j < 12) {
      rslt_soft <- train_test_lm_clust_soft_seq(dist_matrix=dists[[1]], nClusts=j, dat=mvad_covars, Y.cont=num_month_em_last_year, train_idx=train_idx, test_idx=test_idx, fuzziness = 1.5, train_seq_pcs=train_k_scores, test_seq_pcs = test_k_scores)
      rmse.cv.om_seq_soft[i,j] <- rslt_soft$rmse
    }
  }
}


### CFDA 

M <- 36-1 #number of initial months-1

mvad_wide <- mvad[,c(1,15:(15+M))] # Goes only to June 96 
colnames(mvad_wide) <- c("id", as.character(0:M))

# have time and id combinations (so basically only one state column)
mvad_long <- gather(mvad_wide, key = time, value = state, "0":as.character(M))
mvad_long$time <- as.numeric(mvad_long$time)
mvad_long <- mvad_long %>% mutate(state=as.factor(state))
summary_cfd(mvad_long)

# first element is the range overwhich the functions can be evaluated 
# nbasis - number of basis functions (here this is 6 because we have 6 states)
# norder - the order of their degree (here cubic - why cubic?)
basis <- create.bspline.basis(c(0, M), nbasis = 6, norder = 4)




for (i in 1:folds) {  
  cat("Fold Number ",i,"\n")
  
  # pull data
  test_idx <- idx == i  
  train_idx <- !test_idx
  mvad_long_train <- mvad_long %>% filter(id %in% ids[train_idx])
  mvad_long_test <- mvad_long %>% filter(id %in% ids[test_idx])
  
  # compute encodings for - get warning messages that at least one states not in support of 
  # one basis function (I assume this is because of HE in the earlier years)
  fmca.train <- compute_optimal_encoding(mvad_long_train, basis, nCores = 7,verbose=F)
  pcs.train <- fmca.train$pc
  nComps <- ncol(pcs.train)
  
  colnames(pcs.train) <- paste0("PC",1:nComps)
  pcs.test <- predict(fmca.train,newdata=mvad_long_test,method="parallel",nCores=7)
  colnames(pcs.test) <- paste0("PC",1:nComps)
  
  #it's relatively easy to run this with 5 or 10 harmonics - or in the equation, below
  train_harm <- mvad_covars[train_idx,] %>% add_column(as_tibble(pcs.train)) %>% mutate(num_month_em_last_year=num_month_em_last_year[train_idx])
  test_harm <- mvad_covars[test_idx,] %>% add_column(as_tibble(pcs.test)) %>% mutate(num_month_em_last_year=num_month_em_last_year[test_idx])
  
  for (j in 0:nHarms) {
    # create equation
    fmla.str <- paste0("num_month_em_last_year~",paste0(names(mvad_covars),collapse="+"))
    
    # add in pcs to equation, based on the number of cfda pcs to use 
    if (j>0) {
      for (k in 1:j) {
        fmla.str <- paste0(fmla.str,"+PC",k)
      }
    }
    # named formula dictates the model run! 
    fmla <- as.formula(fmla.str)
    mvad_lm_harm <- lm(fmla, data=train_harm)
    pred_harm <- predict(mvad_lm_harm, newdata = test_harm)
    
    rmse.cv.harm[i,j] <- sqrt(sum((pred_harm - test_harm$num_month_em_last_year)^2)/ nrow(test_harm))
  }
} 


## Competition Evaluation
# Ultimately - we want the data and a plot 

om_hard_means_lm <- apply(rmse.cv.om_hard, 2, mean) 
om_soft_means_lm <- apply(rmse.cv.om_soft, 2, mean)
lcs_hard_means_lm <- apply(rmse.cv.lcs_hard, 2, mean)
lcs_soft_means_lm <- apply(rmse.cv.lcs_soft, 2, mean)
hard_cluster_seqs_lm <- apply(rmse.cv.om_seq_hard, 2, mean)
soft_cluster_seqs_lm <- apply(rmse.cv.om_seq_soft, 2, mean)
windows_means_lm <- apply(rmse.cv.windows, 2, mean)
harm_means_lm <- apply(rmse.cv.harm, 2, mean)


rmse.frame <- data.frame(matrix(NA, nrow = 25, ncol = 8))
colnames(rmse.frame) <- c("OM_hard", "OM_soft", "LCS_hard", "LCS_soft", "Windows", "CFDA", "OM_hard_SeqPC", "OM_soft_SeqPC")

rmse.frame$OM_hard <- om_hard_means_lm
rmse.frame$OM_soft <- om_soft_means_lm
rmse.frame$LCS_hard <- lcs_hard_means_lm
rmse.frame$LCS_soft <- lcs_soft_means_lm
rmse.frame$Windows <- c(windows_means_lm, rep(NA,9))
rmse.frame$CFDA <- harm_means_lm
rmse.frame$OM_hard_SeqPC <- hard_cluster_seqs_lm
rmse.frame$OM_soft_SeqPC <- soft_cluster_seqs_lm
rmse.frame$Index <- 1:nrow(rmse.frame)

rmse_long <- rmse.frame %>%
  pivot_longer(
    cols = -Index,          # Columns to pivot (all except Index)
    names_to = "Method",    # New column for the column names (OM, LCS, etc.)
    values_to = "RMSE"      # New column for the values
  )

ggplot(data = rmse_long, 
       aes(x = Index, y = RMSE, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  labs(title = "Linear Regression Performance by Method",
       x = "Components/Clusters",
       y = "RMSE Value",
       color = "Method") +
  theme_minimal()


write.csv(rmse.frame, "LinearCompetitionRMSE.csv", row.names=FALSE)
pdf("LinearCompPlot.pdf",width=8,height=6)


