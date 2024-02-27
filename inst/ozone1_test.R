# Getting the data
# install.packages("SemiPar")
library(earth)
library(tidyverse)
rm(list=ls())
devtools::load_all("/users/research/mmarques/spline_bart_lab/rspBART26/")
set.seed(42)
# Loading the data
type_ <- "ozone1"
data("ozone1")

# Selecting specific continuous variables
real_data <- ozone1 %>% dplyr::select(O3,temp,humidity) %>% dplyr::filter(humidity>min(humidity))
real_data <- real_data %>% dplyr::rename(y = O3)

# Setting the training and a test set
train <- real_data
test_index <- sample(1:nrow(real_data), size = 50)
test <- train[test_index,,]

# Splitting the x.data component and the y.data component
x_train <-train[,colnames(train)!="y", drop = FALSE]
y_train <- train$y
x_test <- test[,colnames(train)!="y",drop = FALSE]


# Setting the main components of the model
nIknots <- 10
dif_order <- 2
n_tree <- 10
n_mcmc <- 5000
n_burn <- 3000
node_min_size <- 25

rsp_mod <- rspBART(x_train = x_train,y_train = y_train,
                   x_test = x_test, nIknots = nIknots,
                   dif_order = dif_order, n_tree = n_tree,
                   n_mcmc = n_mcmc,node_min_size = node_min_size,
                   varimportance_bool = TRUE,
                   main_effects_pred = TRUE,alpha = 0.5)


saveRDS(object = rsp_mod,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART26/",type_,"/single_run/v31_intercept_single_run_nIknots_",
                                       nIknots,"_ntree_",n_tree,"_nodesize_",node_min_size,".Rds"))
