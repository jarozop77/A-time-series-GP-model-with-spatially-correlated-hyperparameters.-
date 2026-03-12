rm(list = ls())

# Libraries -----

library(tidyverse)
library(sf)
library(spdep)
library(ggplot2)

# Inputs ---- 

#Malaria data
malaria_data = readRDS(file = "data/dataset_toyexample.RDS")
summary(malaria_data)

#Shapefile 

shp_file = readRDS(file = "data/shape_file.RDS")

#Starting values
cov.pars = readRDS("data/input_files_init_values/param_GPmle_models_mod.RDS")
cov_cor_indM = readRDS("data/input_files_init_values/cov_mat_GPmle_mod.RDS")

# Data wrangling ----

shp_file = shp_file %>% 
  dplyr::arrange(distcode)
matrix.adj <- poly2nb(shp_file,queen=TRUE)

# Step : find districts with no neighbors
no_nb <- which(card(matrix.adj) == 0)

# Step 3: compute centroids
centroids <- st_coordinates(st_centroid(shp_file))

# Step 4: for each isolated district, assign closest neighbor
for(i in no_nb){
  
  # distances from district i to all others
  dists <- sqrt(
    (centroids[,1] - centroids[i,1])^2 +
      (centroids[,2] - centroids[i,2])^2
  )
  
  dists[i] <- Inf  # ignore itself
  
  # find nearest district
  nearest <- which.min(dists)
  
  # assign as neighbor (both directions)
  matrix.adj[[i]] <- nearest
  matrix.adj[[nearest]] <- unique(c(matrix.adj[[nearest]], i))
}

data = malaria_data %>% 
  dplyr::arrange(distcode)

idx = unique(data$distcode)

cov.pars.filter = cov.pars %>% 
  dplyr::filter(distcode %in% idx) %>% 
  dplyr::arrange(distcode)

cov_cor_indM_sub <- cov_cor_indM[as.character(idx)]

# parameters ----# 

set.seed(7021)
ntraining = 53
ntest = 7 

var_area_id = "distcode"
var_t = "t"
niter = 5000

# sampler ---- 

source("scripts/Inputs_sampler.R")

inputs <- input_user_sampler(
  data = data, var_reg = "distcode", var_t = "t", 
  X = data %>% dplyr::select(sin12,cos12), 
  adj_mat = matrix.adj,
  start_values = list(logsigma = cov.pars.filter$logsigma,
                      logphi = cov.pars.filter$logphi,
                      lognu = cov.pars.filter$lognu), 
  train_test_point = 82, 
  var_y = "log_inc"
)

source("scripts/sampler.R")

timing <- system.time({
  sampler_mcmc <- sampler_mcmc_array(niter)
})

timing

saveRDS(sampler_mcmc, file = "results/Results_toyexample.RDS")

