# Libraries ---------------------------------------------------------------

#For spatial analysis 
#library(spdep)       
library(sf)
library(sp)
library(spatialreg)

#For modelling purposes
library(cmdstanr)
library(posterior)

#For wrangling dataset & visualization
library(tidyverse)
library(ggplot2)

#For the variogram
library(geoR)

#For splines
library(splines)

# Function preparing dataset  ---------------------------------------------

preparing_data_function = function(data, t_test, distribution){
  if(distribution == "GAUSSIAN"){
  data = data %>%
  dplyr::arrange(distcode, year, month) %>%
  dplyr::group_by(distcode) %>%
  dplyr::mutate(cases = ifelse(is.na(cases), 0, cases),
         total_population = ifelse(is.na(total_population), 0, total_population)) %>% 
  dplyr::mutate(t = row_number(), inc = cases/total_population, 
         log_inc = log(inc+0.0001), 
         districtID = cur_group_id()) %>%
  dplyr::ungroup() 
  
  table(data$t)
  data_training = data %>% filter(t<= t_test)
  t_min = min(data_training$t, na.rm = TRUE)
  t_max = max(data_training$t, na.rm = TRUE)
  t_span = t_max - t_min
  y_mean = mean(data_training$log_inc, na.rm = TRUE)
  y_sd = sd(data_training$log_inc, na.rm = TRUE)
  
  data = data %>% 
    dplyr::mutate(uno = 1,
         t_stad = (((2*(t-t_min))/(t_max-t_min))-1),
         log_inc_stad = ((log_inc-y_mean)/y_sd), 
         sin12 = sin(pi*1*t_stad*((t_span)/12)), 
         cos12 = cos(pi*1*t_stad*((t_span)/12)),
         sin6 = sin(pi*2*t_stad*((t_span)/12)), 
         cos6 = cos(pi*2*t_stad*((t_span)/12))) %>% 
    dplyr::arrange(districtID,t) 

  training_data = data %>% dplyr::filter(t <= t_test)
  test_data = data %>% dplyr::filter(t > t_test)
  }
  else {
    data = data %>%
      dplyr::arrange(distcode, year, month) %>%
      dplyr::group_by(distcode) %>%
      dplyr::mutate(cases = ifelse(is.na(cases), 0, cases),
                    total_population = ifelse(is.na(total_population), 0, total_population)) %>% 
      dplyr::mutate(t = row_number(), 
                    districtID = cur_group_id()) %>%
      dplyr::ungroup() 
    
    table(data$t)
    t_min = min(data$t, na.rm = TRUE)
    t_max = max(data$t, na.rm = TRUE)
    t_span = t_max - t_min
    
    data = data %>% 
      dplyr::mutate(uno = 1,
                    t_stad = (((2*(t-t_min))/(t_max-t_min))-1),
                    sin12 = sin(pi*1*t_stad*((t_span)/12)), 
                    cos12 = cos(pi*1*t_stad*((t_span)/12)),
                    sin6 = sin(pi*2*t_stad*((t_span)/12)), 
                    cos6 = cos(pi*2*t_stad*((t_span)/12))) %>% 
      dplyr::arrange(districtID,t) 
    
    training_data = data %>% dplyr::filter(t <= t_test)
    test_data = data %>% dplyr::filter(t > t_test)
  }
  return(list(training_data, test_data))
}


# Function data for Stan --------------------------------------------------

data_4_stan_function = function(data, shp_file, n_vars, n_splines=0, 
                                splines, distribution){

  #Parameters
  ntime = length(unique(data$t))
  ndistricts = length(unique(data$distcode))
  nsplines = n_splines
  
  # Input iCAR
  coords <- sf::st_coordinates(sf::st_centroid(shp_file))
  matrix.adj <- st_touches(shp_file)
  
  # Find isolated regions
  
  cardinality <- lengths(matrix.adj)   # number of neighbors per polygon
  islands <- which(cardinality == 0) 
  
  for (i in islands) {
    # Compute distances to all others
    dist_vec <- sp::spDistsN1(coords, coords[i, ], longlat = FALSE)
    # Find nearest neighbor index (excluding self)
    nearest <- order(dist_vec)[2]
    matrix.adj[[i]] <- nearest
    # Make the connection symmetric
    matrix.adj[[nearest]] <- unique(c(matrix.adj[[nearest]], i))
  }
  
  n <- length(matrix.adj)
  W.matrix <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (length(matrix.adj[[i]]) > 0) W.matrix[i, matrix.adj[[i]]] <- 1
  }
  D <- diag(rowSums(W.matrix))
  
  # Convert adjacency matrix to edge list
  edges <- which(W.matrix == 1, arr.ind = TRUE)
  
  # Keep only upper triangular (avoid double-counting i–j and j–i)
  edges <- edges[edges[,1] < edges[,2], ]
  
  # Rule of thumb grainsize
  #threads <- parallel::detectCores()   # or manually set
  #grainsize <- max(1, floor(ndistricts / (2 * threads)))
  
  # Stan inputs
  n_edges <- nrow(edges)
  node1 <- edges[,1]
  node2 <- edges[,2]
  
  #Input covariates 
  
  sincos_terms_training = vector(mode = "list", length = ndistricts)
  
  for(i in seq_len(ndistricts)){
    sincos_terms_training[[i]] = data %>% 
      dplyr::filter(districtID == i) %>% 
      dplyr::select(sin12,cos12,sin6,cos6) %>% as.matrix()
  }
  
  # Input splines
  time_stan_model = unique(data$t_stad)
  t <- time_stan_model
  tmin <- min(t)
  tmax <- max(t)
  
  if(splines == TRUE){
  # Use equally spaced internal knots
  internal_knots <- seq(tmin, tmax, length.out = n_splines - 2)
  range_splines = seq(tmin,tmax,length.out = n_splines)
  B <- bs(t,
          knots = internal_knots,
          intercept = FALSE,
          Boundary.knots = c(tmin, tmax))[,-(n_splines + 1)]
  
  all_knots <- c(tmin, internal_knots, tmax)
  h <- median(diff(all_knots))
  }
  
  if(distribution == "GAUSSIAN" & splines == TRUE){
    #Input y and time
    log_inc = data$log_inc_stad
    data_log_inc <- t(matrix(log_inc, nrow = ntime, ncol = ndistricts))
    
    data_districts_stan <- list(n_times = ntime, n_sites = ndistricts,
                                n_sincos = n_vars, n_edges = n_edges,
                                n_splines = nsplines, spline_time = B,
                                splines_range = range_splines,
                                y = data_log_inc, 
                                seasonality = sincos_terms_training, 
                                node1 = node1, node2=node2)
  }
  
  if(distribution == "GAUSSIAN" & splines == FALSE){
    #Input y and time
    log_inc = data$log_inc_stad
    data_log_inc <- t(matrix(log_inc, nrow = ntime, ncol = ndistricts))
    
    data_districts_stan <- list(n_times = ntime, n_sites = ndistricts,
                                n_sincos = n_vars, n_edges = n_edges,
                                y = data_log_inc, time = time_stan_model,
                                seasonality = sincos_terms_training, 
                                node1 = node1, node2=node2)
  }
  
  if(distribution %in% c("POISSON","BINOMIAL") & splines == FALSE){
    #Input y and time
    cases = data$cases
    data_cases <- t(matrix(cases, nrow = ntime, ncol = ndistricts))
    
    population_training = vector(mode = "list", length = ndistricts)
    for(i in seq_len(ndistricts)){
      pop = data %>% 
        dplyr::filter(districtID == i) %>% 
        dplyr::select(total_population) 
      population_training[[i]] = pop$total_population
    }
    
    data_districts_stan <- list(n_times = ntime, n_sites = ndistricts,
                                n_sincos = n_vars, n_edges = n_edges,
                                y = data_cases, time = time_stan_model,
                                seasonality = sincos_terms_training,
                                population = population_training, 
                                node1 = node1, node2=node2)
  }
  if(distribution %in% c("POISSON","BINOMIAL") & splines == TRUE){
    #Input y and time
    cases = data$cases
    data_cases <- t(matrix(cases, nrow = ntime, ncol = ndistricts))
    
    population_training = vector(mode = "list", length = ndistricts)
    for(i in seq_len(ndistricts)){
      pop = data %>% 
        dplyr::filter(districtID == i) %>% 
        dplyr::select(total_population) 
      population_training[[i]] = pop$total_population
    }
    
    data_districts_stan <- list(n_times = ntime, n_sites = ndistricts,
                                n_sincos = n_vars, n_edges = n_edges,
                                y = data_cases, time = time_stan_model,
                                n_splines = nsplines, spline_time = B,
                                splines_range = range_splines,
                                seasonality = sincos_terms_training,
                                population = population_training, 
                                node1 = node1, node2=node2)
  }
  
  
  return(data_districts_stan)
}

#dt_stan = data_4_stan_function(dt[[1]], sho_data)
  



