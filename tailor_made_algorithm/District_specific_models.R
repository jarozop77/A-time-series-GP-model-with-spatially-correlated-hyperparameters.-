rm(list = ls())

# Libraries ---------------------------------------------------------------

library(tidyverse) #To manipulate data. 
library(sf) #to plot maps
library(stringi) #To transform especial characters
library(zoo)
library(geoR)
library(PrevMap)

# Datasets ----------------------------------------------------------------

malaria_data = readRDS("data/dataset_toyexample.RDS")

# Parameters --------------------------------------------------------------

ndistricts = length(unique(malaria_data$distcode))
kp = 1.5
ntime = 53
t_training = 53

set.seed(2107)

training_data = malaria_data %>% filter(t<=t_training)
test_data = malaria_data %>% filter(t>t_training)

# Filtering the datasets --------------------------------------------------

dIDs_table = training_data %>% 
  dplyr::select(distcode) %>% unique() %>% 
  dplyr::arrange(distcode) %>% 
  dplyr::mutate(districtID = row_number())
training_data = training_data %>% left_join(dIDs_table)

distri_datasets = vector(mode = "list", length = ndistricts)

list_data = function(data, dID){
  distr_data = data %>%
    dplyr::filter(districtID == dID) %>% 
    dplyr::mutate(uno = 1)
  return(distr_data)
}

for(i in 1:ndistricts){
  distri_datasets[[i]] = list_data(training_data, i)
}

# Plot of the incidence or counts -----------------------------------------

distri_descrplot = vector(mode = "list", length = ndistricts)

descriptive_plot = function(data, var){
  
  title = paste0(data$distname[1], " (", data$districtID[1], ")")
  
  ggplot(data = data, aes(x = date, y = {{ var }})) +
    geom_line(colour = "darkblue") +
    theme(axis.text.x = element_text(size = 8), 
          axis.line = element_line(linewidth = 0.5, colour = "darkblue", linetype = 1)) + 
    ylab("Number of malaria cases (log scale)") + 
    xlab("Year-month") +
    ggtitle(title)
}

for(i in 1:ndistricts){
  distri_descrplot[[i]] = descriptive_plot(distri_datasets[[i]], log_inc)
}

distri_descrplot[[16]]

# Variogram plot ----------------------------------------------------------

variog = function(dataset){
  #dataset = distri_datasets[[7]]
  fit1 = lm(log_inc ~ sin12 + cos12,
            data = dataset)     
  
  res = data.frame(t=as.numeric(names(residuals(fit1))),
                   residuals1 = residuals(fit1))
  
  dataset = dataset %>% 
    left_join(res)
  
  #Step 2.2: Variogram
  data = dplyr::filter(dataset, is.na(residuals1)==FALSE)
  coordinates(data) = ~t + uno
  
  variogram1 = geoR::variog(data,
                            coords = data@coords, 
                            data = data$residuals1,
                            uvec = seq(1,45,1))
  env1 = variog.mc.env(data, coords = data@coords, data = data$residuals1,
                       obj.variog = variogram1, nsim = 100)
  info_vario1 = data.frame(cbind(u=env1$u,vl=env1$v.lower,v=variogram1$v,vu=env1$v.upper))
  
  title = paste0("Variogram ", unique(dataset$dist_name), " (", unique(dataset$districtID), ")")
  
  plot = ggplot(info_vario1, aes(x=u, y=v)) + geom_line() + 
    geom_line(aes(y=vl), color="steelblue") +
    geom_line(aes(y=vu), color="steelblue") + ylab("Variogram") + 
    xlab("Time separation (months)") + ggtitle(title)
  
  return(plot)
  
}

variogram_plot = vector(mode = "list", length = ndistricts)

for(i in 1:ndistricts){
  variogram_plot[[i]] = variog(distri_datasets[[i]])
}
variogram_plot[[16]]

# Data for modelling without NA -------------------------------------------

data.formodel = function(dataset){
  #dataset = distri_datasets[[1]]
  fit1 = lm(log_inc ~ sin12 + cos12,
            data = dataset)      
  
  res = data.frame(t=as.numeric(names(residuals(fit1))),
                   residuals1 = residuals(fit1))
  
  dataset = dataset %>% left_join(res)
  
  #Step 2.2: Variogram
  data = dplyr::filter(dataset, is.na(residuals1)==FALSE)
  
  return(data)
}

d4_model = vector(mode = "list", length = ndistricts)

for(i in 1:ndistricts){
  d4_model[[i]] = data.formodel(distri_datasets[[i]])
}

# Initial values  ---------------------------------------------------------

ini_vals15 = vector(mode = "list", length = ndistricts)

initial_values = function(dataset, kappa){
  coordinates(dataset) = ~t + uno
  #Variogram
  variogram1 = geoR::variog(dataset, coords = dataset@coords, 
                            data = dataset$residuals1,
                            uvec = seq(1,24,1))
  
  #Initial values for the variofit
  sigma.vals = seq(min(variogram1$v), max(variogram1$v), by=0.01)
  nug.vals = seq(min(variogram1$v), max(variogram1$v), by=0.01)
  phi.vals = seq(1,5,by = 1)
  
  ini.vals = expand.grid(sigma.vals, phi.vals)
  
  vari.fit <- variofit(variogram1, ini.cov.pars = ini.vals,
                       cov.model = "matern",
                       fix.nugget = FALSE, nugget = nug.vals,
                       fix.kappa = TRUE, kappa = kappa)
  #Validating the values 
  tau=ifelse(vari.fit$nugget==0, min(variogram1$v), vari.fit$nugget)
  phi=ifelse(vari.fit$cov.pars[2]==0, 1, vari.fit$cov.pars[2])
  sigma=ifelse(vari.fit$cov.pars[1]==0, min(variogram1$v), vari.fit$cov.pars[1])
  #sigma=ifelse(vari.fit$cov.pars[1]>max(variogram1$v), max(variogram1$v), vari.fit$cov.pars[1])
  
  values = c(tau, sigma, phi)
  return(values)
}

for(i in 1:ndistricts){
  ini_vals15[[i]] = initial_values(d4_model[[i]], kappa=kp)
}

ini_vals15[[1]] 

# Model -------------------------------------------------------------------

model = function(dataset, val_ini, kappa){
  fit <- linear.model.MLE(form= log_inc ~ sin12 + cos12,
                          coords = ~t + uno,
                          start.cov.pars = c(val_ini[3],
                                             val_ini[1]/val_ini[2]),
                          kappa=kappa,
                          data=dataset,
                          method="nlminb")
  return(fit)
}

fit_k15 = vector(mode = "list", length = ndistricts)

for(i in 1:ndistricts){
  fit_k15[[i]] = model(d4_model[[i]], ini_vals15[[i]], kappa=kp)
  names(fit_k15)[i] = as.character(d4_model[[i]] %>% dplyr::select(distcode) %>% distinct(distcode) %>% pull())
}

summary(fit_k15[[7]])

# Model-based estimates  --------------------------------------------------

estimates = function(model, dataset){
  # i=7
  # model = fit_k15[[i]]
  # dataset = distri_datasets[[i]]
  
  predx = as.data.frame(dataset %>% dplyr::select(uno, sin12, cos12))
  
  pred_or <- spatial.pred.linear.MLE(object = model,
                                     predictors = predx,
                                     grid.pred = cbind(seq(1,ntime,1), rep(1, ntime)),
                                     scale.predictions = "odds")
  
  pred_log <- spatial.pred.linear.MLE(object = model,
                                      predictors = predx,
                                      grid.pred = cbind(seq(1,ntime,1), rep(1, ntime)),
                                      scale.predictions = "logit")
  
  dataset$est.m4 = pred_log$logit$predictions
  dataset$est.m4.ll = pred_log$logit$quantiles[,1]
  dataset$est.m4.ul = pred_log$logit$quantiles[,2]
  
  return(dataset)
}

data_est = vector(mode = "list", length = ndistricts)

for(i in 1:ndistricts){
  data_est[[i]] =  try(estimates(fit_k15[[i]],distri_datasets[[i]]))
}

data_est_res = data_est %>% bind_rows()

# Predictions -------------------------------------------------------------

test_data = test_data %>% 
  dplyr::arrange(t)

distcodes = unique(training_data$distcode)

#Dataset for predictions 
test_datasets_ds = vector(mode = "list", length = ndistricts)

for(i in 1:ndistricts){
  test_datasets_ds[[i]] = test_data %>% 
    dplyr::filter(distcode==distcodes[i]) %>% 
    dplyr::mutate(uno = 1)
}

predictions = function(model, data){
  # model = fit_k15[[1]]
  # data = test_datasets_ds[[1]]
  ntest = length(unique(data$t))
  
  predx = as.data.frame(data %>% ungroup() %>% dplyr::select(uno, sin12, cos12))
  grid = as.matrix(cbind(data$t, rep(1, ntest)))
  pred0.5 <- spatial.pred.linear.MLE(object = model,
                                     predictors = predx,
                                     grid.pred = grid,
                                     scale.predictions = "logit")
  
  data$est.m4 = pred0.5$logit$predictions
  data$est.m4.ll = pred0.5$logit$quantiles[,1]
  data$est.m4.ul = pred0.5$logit$quantiles[,2]
  #data$est.var.m2 = abs(data$pred.m2 - data$pred.m2.ll)/1.959
  return(data)
}

distri_pred = vector(mode = "list", length = ndistricts)

for(i in 1:ndistricts){
  distri_pred[[i]] = predictions(fit_k15[[i]],test_datasets_ds[[i]])
}

data_pred = distri_pred %>% bind_rows()

# Gathering the fixed and covariance parameters ---------------------------

cov.pars = data.frame(distcode = integer(ndistricts), dist_name = NA,
                      logsigma = NA, se.logsigma = NA, 
                      logphi = NA, se.logphi = NA, 
                      logtau = NA, se.logtau = NA,
                      lognu = NA, se.lognu = NA)

for(i in 1:ndistricts){
  cov.pars$distcode[i] = unique(distri_pred[[i]]$distcode)
  cov.pars$dist_name[i] = unique(distri_pred[[i]]$distname)
  cov.pars$logsigma[i] = summary(fit_k15[[i]])$cov.pars[1,1]
  cov.pars$se.logsigma[i] = summary(fit_k15[[i]])$cov.pars[1,2]
  cov.pars$logphi[i] = summary(fit_k15[[i]])$cov.pars[2,1]
  cov.pars$se.logphi[i] = summary(fit_k15[[i]])$cov.pars[2,2]
  cov.pars$logtau[i] = summary(fit_k15[[i]])$cov.pars[3,1]
  cov.pars$se.logtau[i] = summary(fit_k15[[i]])$cov.pars[3,2]
  cov.pars$lognu[i] = cov.pars$logtau[i] - cov.pars$logsigma[i]
}

cov_cor_indM = vector(mode = "list", length = ndistricts)

for(i in 1:ndistricts){
  mat = fit_k15[[i]]$covariance[4:6,4:6]
  cov_cor_indM[[i]] = mat
  cov.pars$se.lognu[i] = cov_cor_indM[[i]][3,3]
  names(cov_cor_indM)[i] = names(fit_k15)[i]
}

cov_cor_indM[[15]]

# cov_cor_indM: named list of 3x3 matrices

check_and_fix_vcov <- function(lst, lower = -4, upper = 4) {
  changed <- logical(length(lst))
  
  fixed <- lapply(seq_along(lst), function(i) {
    M <- lst[[i]]
    
    # flag if any element is NA/NaN/Inf or outside bounds
    bad <- any(!is.finite(M)) || any(M < lower | M > upper)
    
    changed[i] <<- bad
    if (bad) diag(nrow(M)) else M
  })
  
  # keep names
  names(fixed) <- names(lst)
  
  list(
    fixed_list = fixed,
    n_changed = sum(changed),
    changed_index = which(changed),
    changed_names = names(lst)[changed]
  )
}
res <- check_and_fix_vcov(cov_cor_indM, lower = -4, upper = 4)

cov_cor_indM_fixed <- res$fixed_list
res$n_changed          # how many matrices were replaced
res$changed_names      # which ones (optional)

# helper: replace NA/NaN/Inf, and replace out-of-bounds values with a default
fix_bounds <- function(x, lower, upper, default) {
  x[!is.finite(x)] <- default
  x[x < lower | x > upper] <- default
  x
}

# helper: cap SEs above a threshold (and also fix non-finite)
fix_se <- function(x, max_se, default) {
  x[!is.finite(x)] <- default
  x[x > max_se] <- default
  x
}

cov.pars <- cov.pars %>%
  mutate(
    across(c(lognu, logsigma, logphi, logtau), ~ fix_bounds(.x, lower = -6, upper = 4, default = -2)),
    across(c(se.lognu, se.logsigma, se.logphi, se.logtau), ~ fix_se(.x, max_se = 4, default = 0.21))
  )

saveRDS(cov_cor_indM, file = "data/input_files_init_values/cov_mat_GPmle_mod.RDS")
saveRDS(cov.pars, file = "data/input_files_init_values/param_GPmle_models_mod.RDS")



