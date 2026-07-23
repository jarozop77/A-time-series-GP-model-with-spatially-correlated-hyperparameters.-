setwd("/Users/alejandrorozo/Library/CloudStorage/GoogleDrive-alejandro.rozoposada@uhasselt.be/My Drive/Malaria Project/Malaria/MALE/GitHub_repository/reproducible_code")

rm(list=ls())

# libraries ---- 

library(fs)

# input files ---- 

# Code of the model implementation in stan 
stan_file = "stan_implementation/model/model.stan"

# Script containing preprocessing functions for the data
source("stan_implementation/Preprocessing.R")

# data
data = readRDS("data/simulated_data.RDS")
data = data %>% 
  dplyr::rename(cases = sim.cases)

# shapefile 

shp = readRDS("data/district_boundaries.RDS")

# Parameters for STAN model ----
n_chains <- 1
n_iter <- 100
n_burnin <- 100
n_threads_per_chain <- 1

# Compile STAN model ------------------------------------------------------

mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# Obtaining inputs for STAN model --------------------------------------------

# 1) Training and test datasets
t_test <- 55
prep_out <- preparing_data_function(data, t_test, "GAUSSIAN")

training_data <- prep_out[[1]]
test_data <- prep_out[[2]]

ggplot(shp) + geom_sf()

# 2) DATA 4 STAN 
data_stan <- data_4_stan_function(training_data, shp, 
                                  n_vars = 4, n_splines = 0, 
                                  splines = FALSE,distribution = "GAUSSIAN")
data_stan$grainsize <- 1

# Compilation of the model  -----------------------------------------------

fit <- mod$sample(
  data = data_stan,
  #init = list(init_list),
  chains = n_chains,
  iter_sampling = n_iter,
  iter_warmup = n_burnin,
  parallel_chains = 1,         
  threads_per_chain = 10,
  refresh = 5, #How often it is printed
  adapt_delta = 0.95,
  max_treedepth = 15,
  output_dir = "stan_implementation"
)

# Save outputs ----------------------------------

saveRDS(training_data, file = "outputs/stan_implementation/training_data.RDS")
saveRDS(test_data, file = "outputs/stan_implementation/test_data.RDS")

# Save cmdstan fit object 
fit$save_object(file = "outputs/stan_implementation/model_fitted.RDS")




