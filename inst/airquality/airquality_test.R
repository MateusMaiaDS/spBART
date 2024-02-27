# Loading the packages and everything else needed to
library(tidyverse)
rm(list=ls())

# This file is to compare first trial of predictive performance of spBART and BART
library(dbarts)
library(mlbench)
library(purrr)
library(MOTRbart)
library(doParallel)
devtools::load_all("/users/research/mmarques/spline_bart_lab/rspBART26/")
source("/users/research/mmarques/spline_bart_lab/rspBART26/R/sim_functions_26.R")
source("/users/research/mmarques/spline_bart_lab/rspBART26/R/main_function_26.R")

# Setting the seed
seed_ <- 42
set.seed(seed_)
competitors_only <- FALSE

# Setting the data
airquality <- read_csv("/users/research/mmarques/spline_bart_lab/rspBART26/inst/airquality/airquality.csv")
# Removing the NA columns
airquality <- airquality[complete.cases(airquality),]

# Transforming the data so it is aligned with the Bayesian MARS paper
data_ <- airquality %>% dplyr::select(Ozone,Solar.R,Wind,Temp) %>% dplyr::mutate(Ozone = (Ozone)^(1/3)) %>% dplyr::rename(y = Ozone)


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
node_min_size_ <- 25
n_mcmc_ <- 10000
n_burn_ <- 5000
pen_basis_ <- TRUE

# Selecting a simulated scenarion
type_ <- c("airquality")

# ================
# Printing message
# ================

print(paste0("Airquality dataset => nIknots: ", nIknots_,
             " Ntree: ",ntree_, " Seed: ",seed_, " Alpha:", alpha_,
             "Update \tau_\beta: ", update_tau_beta_, "_type_", type_,
             "_nmcmc_",n_mcmc_, "_nburn_",n_burn_))


cv_ <- vector("list", n_rep_)

# Getting the cv element
cv_<- kfold(data_ = data_,nfold_ = n_rep_,seed_ = seed_)

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
    aux <- competitors_comparison_(cv_element= cv_[[i]],
                                   fold_ = i,seed_ = seed_,
                                   return_models = FALSE)
  }
  # }

  aux
}


stopCluster(cl)




# Saving the plots
if(competitors_only){
  saveRDS(object = result,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART26/",type_,"/competitors_seed_",seed_,".Rds"))
} else {
  saveRDS(object = result,file = paste0("/users/research/mmarques/spline_bart_lab/preliminar_results/rspBART26/",type_,"/v31_intercept_psBART_seed_",seed_,
                                        "_nIknots_",nIknots_,"_ntree_",ntree_,
                                        "_alpha_",alpha_,"_dif_",dif_order_,"_nmin_",node_min_size_,
                                        "_nmcmc_",n_mcmc_,"_nburn_",n_burn_,".Rds"))
}




