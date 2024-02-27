# This file is to compare first trial of predictive performance of spBART and BART
rm(list=ls())
library(dbarts)
library(mlbench)
library(purrr)
library(MOTRbart)
library(doParallel)
source("/users/research/mmarques/spline_bart_lab/rspBART26/R/sim_functions_26.R")
source("/users/research/mmarques/spline_bart_lab/rspBART26/R/main_function_26.R")
set.seed(42)
competitors_only <- TRUE


n_ <- 1000
sd_ <- 1
n_rep_ <- 10
nIknots_ <- 20
ntree_ <- 10
dif_order_ <- 2
use_bs_ <- FALSE
seed_ <- 42
y_scale_ <- TRUE
motr_bart_ <- FALSE
all_ <- FALSE
alpha_ <- 0.5
stump_ <- FALSE
scale_init_ <- FALSE
update_tau_beta_ <- TRUE
node_min_size_ <- 50
n_mcmc_ <- 5000
n_burn_ <- 3000
pen_basis_ <- TRUE

# Selecting a simulated scenarion
# (1): "oned_break" one dimensionnal sin(2*x) with a break
# (2): "friedman_nointer_nonoise": four-dimensional friedmna setting with no interaction terms and no extra X noise variables
# (3): "interaction
type_ <- c("friedman")
# type_ <- c("friedman_break")
# type_ <- "smooth.main.formula"
# type_ <- "non.smooth.main.formula"
# type_ <- "non.and.smooth.main.formula"
# type_ <- 'mlbench.d1.break'

# ================
# Printing message
# ================

print(paste0("N: ",n_," SD: ", sd_, " nIknots: ", nIknots_,
             " Ntree: ",ntree_, " Seed: ",seed_, " Alpha:", alpha_,
             "Update \tau_\beta: ", update_tau_beta_, "_type_", type_,
             "_nmcmc_",n_mcmc_, "_nburn_",n_burn_))


cv_ <- vector("list", n_rep_)


# Generating CV_ object
for( i in 1:n_rep_){


  # train <- mlbench.d1.break(n = n_,sd = sd_) %>% as.data.frame()
  # test <- mlbench.d1.break(n = n_,sd = sd_) %>% as.data.frame()

  if(type_ == "smooth.main.formula"){
    train_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = smooth.main.formula)
    test_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = smooth.main.formula)
    train <- data.frame(x = train_sim$x, y = train_sim$y)
    test <- data.frame(x = test_sim$x, y = test_sim$y)
  }

  if(type_ == "non.smooth.main.formula"){
    train_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.smooth.main.formula)
    test_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.smooth.main.formula)
    train <- data.frame(x = train_sim$x, y = train_sim$y)
    test <- data.frame(x = test_sim$x, y = test_sim$y)
  }

  if(type_ == "non.and.smooth.main.formula"){
    train_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.and.smooth.main.formula)
    test_sim <- sim.gen(n_ = n_,sd_ = sd_,p_ = 10,formula_ = non.and.smooth.main.formula)
    train <- data.frame(x = train_sim$x, y = train_sim$y)
    test <- data.frame(x = test_sim$x, y = test_sim$y)
  }

  if(type_ == "friedman_nointer_nonoise"){
    train <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.nointeraction(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_break"){
    train <- break.mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
    test <- break.mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_nointer_noise"){
    train <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.nointeraction.noise(n = n_,sd = sd_) %>% as.data.frame()
  }
  # train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]
  # test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]

  if(type_ == "friedman"){
    train <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == "friedman_interaction"){
    train <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()
    test <- mlbench.friedman1.interaction.only(n = n_,sd = sd_) %>% as.data.frame()
  }

  if(type_ == 'mlbench.d1.break'){
    train <- mlbench.d1.break(n = n_,sd = sd_)  |> as.data.frame()
    test <- mlbench.d1.break(n = n_,sd = sd_) |> as.data.frame()
  }
  cv_[[i]]$train <- train
  cv_[[i]]$test <- test
}

# Setting up the parallel simulation
number_cores <- n_rep_
cl <- parallel::makeCluster(number_cores)
doParallel::registerDoParallel(cl)


# Testing the simple n_tree
result <- foreach(i = 1:n_rep_, .packages = c("dbarts","SoftBart","MOTRbart","dplyr")) %dopar%{

  devtools::load_all("/users/research/mmarques/spline_bart_lab/rspBART26/")
  source("/users/research/mmarques/spline_bart_lab/rspBART26/R/sim_functions_26.R")
  source("/users/research/mmarques/spline_bart_lab/rspBART26/R/main_function_26.R")
  source("/users/research/mmarques/spline_bart_lab/rspBART26/R/cv_functions.R")

  if(isFALSE(competitors_only)){
    aux <- all_spbart_lite_interaction(cv_element = cv_[[i]],
                         nIknots_ = nIknots_,ntree_ = ntree_,
                         node_min_size_ = node_min_size_,
                         seed_ = seed_,
                         alpha_ = alpha_,
                         j = i,dif_order_ = dif_order_,
                         y_scale_ = y_scale_,
                         n_mcmc_ = n_mcmc_,n_burn_ = n_burn_,
                         pen_basis_ = pen_basis_)
  } else {
    aux <- competitors_comparison_(cv_object_fold_ = cv_[[i]],
                                   fold_ = i,seed_ = seed_,
                                   return_models = FALSE)
  }
  # }

  aux
}


stopCluster(cl)




# Saving the plots
if(competitors_only){
  saveRDS(object = result,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART26/",type_,"/competitors_n_",n_,
                                        "_sd_",sd_,".Rds"))
} else {
  saveRDS(object = result,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART26/",type_,"/v31_intercept_psBART_n_",n_,
                                        "_sd_",sd_,"_nIknots_",nIknots_,"_ntree_",ntree_,
                                        "_alpha_",alpha_,"_dif_",dif_order_,"_nmin_",node_min_size_,
                                        "_nmcmc_",n_mcmc_,"_nburn_",n_burn_,".Rds"))
}




