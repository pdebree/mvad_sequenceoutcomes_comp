# Comparing Demographics, Clusters and the combination of both. 

library(tidyverse)
library(TraMineR)
library(cluster)
library(cfda)
source("seqout_utils.R")

set.seed(11)

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

# Looking at only the first fold 
test_idx <- idx == 1
train_idx <- !test_idx


train_demos <- mvad_covars[train_idx,] |> mutate(y = num_month_em_last_year[train_idx])
test_demos <- mvad_covars[test_idx,] |> mutate(y = num_month_em_last_year[test_idx])


fit.lin_dem <- lm(y~., data=train_demos) 
summary(fit.lin_dem)
train_demos_lab <- train_demos |> mutate(preds = predict(fit.lin_dem))
test_demos_lab <- test_demos |> mutate(preds = predict(fit.lin_dem, newdata=test_demos))

# create clustering - om-trate with 
clusterward_hard <- agnes(dists[[1]][train_idx,train_idx], diss=TRUE, method="ward")

# fit with 12 cluster (this was optimal in the linear competition)
hard_cluster_data <- hard_cluster(
  clusterward=clusterward_hard, nClusts=12, covars=mvad_covars, 
  y = num_month_em_last_year, train_idx=train_idx, 
  test_idx=test_idx, dist_matrix = dists[[1]])


train_cl <- hard_cluster_data$train_data[,12:13]
test_cl <- hard_cluster_data$test_data[,12:13]

fit.lin_cl <- lm(y~., data=train_cl) 
summary(fit.lin_cl)

train_labs_cl <- train_cl |> mutate(preds = predict(fit.lin_cl))
test_labs_cl <- test_cl |> mutate(preds = predict(fit.lin_cl, newdata=test_cl))

train_sum_cl <- train_labs_cl |> group_by(cluster) |> summarize(mean_preds = mean(preds), var_preds = sd(preds)^2, mean_y = mean(y), var_y = sd(y)^2)
test_sum_cl <- test_labs_cl |> group_by(cluster) |> summarize(mean_preds = mean(preds), var_preds = sd(preds)^2, mean_y = mean(y), var_y = sd(y)^2)


# both together
fit.lin_both <- lm(y~., data=hard_cluster_data$train_data) 
summary(fit.lin_both)

train_labs_both <- hard_cluster_data$train_data |> mutate(preds = predict(fit.lin_both))
test_labs_both <- hard_cluster_data$test_data |> mutate(preds = predict(fit.lin_both, newdata=hard_cluster_data$test_data))

train_sum_both <- train_labs_both |> group_by(cluster) |> summarize(mean_preds = mean(preds), var_preds = sd(preds)^2, mean_y = mean(y), var_y = sd(y)^2)
test_sum_both <- test_labs_both |> group_by(cluster) |> summarize(mean_preds = mean(preds), var_preds = sd(preds)^2, mean_y = mean(y), var_y = sd(y)^2)
