library(invgamma)

### Translation of the dataset to the sampler ###

input_user_sampler = function(data, var_reg, var_t, train_test_point, X, adj_mat, start_values, var_y){
  # validation
  stopifnot(is.data.frame(data))
  stopifnot(is.character(var_reg), length(var_reg) == 1)
  stopifnot(is.character(var_t), length(var_t) == 1)
  stopifnot(is.list(start_values))
  
  # starting values 
  training_data = data %>% dplyr::filter(t<=train_test_point)
  test_data = data %>% dplyr::filter(t>train_test_point)
  
  #safe selection
  area_codes = data %>% dplyr::ungroup() %>% 
    dplyr::select(dplyr::all_of(var_reg)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull()
  
  t = training_data %>% dplyr:: ungroup() %>% 
    dplyr::select(dplyr::all_of(var_t)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull()
  
  if (!is.numeric(t)) stop("`", var_t, "` must be numeric.", call. = FALSE)
  if (!is.numeric(area_codes)) stop("`", var_reg, "` must be numeric.", call. = FALSE)
  
  ntimes = length(t)
  ndistricts = length(area_codes)
  
  #Covariates 
  
  stopifnot(dim(X)[1] == dim(data)[1])
  if(is.null(colnames(X)) || anyNA(colnames(X))){
      var_names = paste0("V",seq_len(ncol(X)))}
  else{var_names = colnames(X)}
  
  X_training = training_data %>% 
    dplyr::select(distcode, all_of(var_names))
  
  X_test = test_data %>% 
    dplyr::select(distcode, all_of(var_names))

  #Adjacency matrix
  W <- nb2mat(adj_mat, style = "B")
  D <- diag(rowSums(W))
  
  Q <- (D - alpha * W)
  
  # ---- validate start_values ----
  
  if (length(start_values) < 3) {
    stop("`start_values` must be a list with at least 3 elements: log.sigma2, log.phi, log.nu2.",
         call. = FALSE)
  }
  
  sv <- start_values[1:3]
  
  # each must be numeric and length == ndistricts
  bad_type <- !vapply(sv, is.numeric, logical(1))
  if (any(bad_type)) {
    stop("The first 3 elements of `start_values` must be numeric vectors.",
         call. = FALSE)
  }
  
  bad_len <- vapply(sv, length, integer(1)) != ndistricts
  if (any(bad_len)) {
    stop("Each of the first 3 elements of `start_values` must have length == ndistricts (",
         ndistricts, "). Got lengths: ",
         paste(vapply(sv, length, integer(1)), collapse = ", "),
         call. = FALSE)
  }
  
  start_df <- dplyr::tibble(
    area_codes  = area_codes,log.sigma2  = sv[[1]],
    log.phi     = sv[[2]],log.nu2     = sv[[3]] )
  
  nbetas = ncol(X) + 1 
  nparameters = ntimes + nbetas + 3
  
  return(list(training_data = training_data, test_data = test_data, 
              area_codes = area_codes, t = t, Q = Q, X_training = X_training, 
              X_test = X_test, start_values = start_df, var_names = var_names,
    ntimes = ntimes, ndistricts = ndistricts,
    nbetas = nbetas, nparameters = nparameters,
    var_y = var_y))
}

#Hyperparameters of the model
delta = rgamma(1, 1, 0.001) #For the prior of the log sigmas
v_phi = rinvgamma(1, 0.1, 0.1) # For the prior of the log phi
lambda_betas = rgamma(1, 0.5, 0.5) #For the priors of the betas
alpha = 0.9 #Estimate for the BYM model

#Parameters sampler
epsilon_0 = 0.01
df_nu = 30
alpha_target_nu = 0.44




