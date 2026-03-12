rm(list = ls())

# --- libraries ----

library(tidyverse)
library(coda)
library(forecast)   
library(bayesplot)
library(patchwork)

# --- inputs ----

#Malaria data
malaria_data = readRDS(file = "data/dataset_toyexample.RDS")
data = malaria_data %>% arrange(distcode, t)
summary(data)

niter = 5000
burnin <- 1000
thin   <- 5
ndistricts = 16
t_training = 53

#Shapefile 

shp_file = readRDS(file = "data/shape_file.RDS")
shp_file = shp_file %>% arrange(distcode)

#mcmc results#mcmc results
mcmc_results = readRDS(file="results/Results_toyexample.RDS")

eps_sigma = mcmc_results$eps_sigma
plot(seq(1,niter,1), eps_sigma)
eps_phi = mcmc_results$eps_phi
plot(seq(1,niter,1), eps_phi[,16])

# ---- data wrangling ----
training_data = data %>% filter(t<=53)
test_data = data %>% filter(t>53)

distcodes = training_data %>% ungroup() %>% 
  dplyr::select(distname, distcode) %>% 
  dplyr::arrange(distcode) %>% 
  unique()

y_list <- lapply(distcodes$distcode, function(id){
  training_data$log_inc[training_data$distcode == id]
})
names(y_list) <- as.character(distcodes$distcode)

X_list <- lapply(distcodes$distcode, function(id){
 training_data[training_data$distcode == id,c("sin12","cos12")]
})
names(X_list) <- as.character(distcodes$distcode)

X_test <- lapply(distcodes$distcode, function(id){
  test_data[test_data$distcode == id,c("sin12","cos12")]
})
names(X_test) <- as.character(distcodes$distcode)

t_training = unique(training_data$t)
t_test = unique(test_data$t)

# --- chains wrangling ----

P <- mcmc_results$P  
S <- mcmc_results$S   

iter_idx <- seq(burnin + 1, dim(P)[1]-1, by = thin)

P_thin <- P[iter_idx, , , drop = FALSE]
S_thin <- S[iter_idx, , , drop = FALSE]   

extract_param <- function(P_arr, districtID, param_name) {
  tibble(iter  = seq_len(dim(P_arr)[1]), 
         value = as.numeric(P_arr[, districtID, param_name]))
}
extract_Sk <- function(S_arr, districtID, k) {
  tibble(iter  = seq_len(dim(S_arr)[1]),
         value = as.numeric(S_arr[, districtID, k]))
}

# --- Geweke diagnostic ----
geweke_trunc_data <- function(x, max_frac = 0.5, step_frac = 0.05,
                              frac1 = 0.1, frac2 = 0.5) {
  x <- as.numeric(x)
  n <- length(x)
  
  max_discard <- floor(max_frac * n)
  step_size   <- max(1L, floor(step_frac * n))
  discards    <- seq(0L, max_discard, by = step_size)
  
  z_scores <- vapply(discards, function(d) {
    xt <- x[(d + 1L):n]
    # if too short, return NA
    if (length(xt) < 20L) return(NA_real_)
    as.numeric(coda::geweke.diag(coda::as.mcmc(xt),
                                 frac1 = frac1, frac2 = frac2)$z)
  }, numeric(1))
  
  tibble(
    discarded = discards,
    discarded_frac = discards / n,
    z = z_scores
  )
}

geweke_trunc_plot <- function(x, param_label = "param",
                              max_frac = 0.5, step_frac = 0.05,
                              frac1 = 0.1, frac2 = 0.5) {
  gdf <- geweke_trunc_data(x, max_frac, step_frac, frac1, frac2)
  
  ggplot(gdf, aes(discarded, z)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed") +
    labs(
      title = paste("Geweke (truncation):", param_label),
      x = "Draws discarded from start",
      y = "Geweke Z"
    ) +
    theme_minimal()
}

# --- Diagnostic plot functions ----
diag_plots_param <- function(P_arr, districtID, param_name, max_lag = 50) {
  # P_arr = P_thin
  # districtID = "900959" 
  # param_name = "logphi" 
  # max_lag = 50
  
  df <- extract_param(P_arr, districtID, param_name)
  
  # trace
  p_trace <- ggplot(df, aes(iter, value)) +
    geom_line() +
    labs(title = paste("Trace:", param_name, "district", districtID),
         x = "Saved draw", y = param_name) +
    theme_minimal()
  
  # histogram
  p_hist <- ggplot(df, aes(value)) +
    geom_histogram(bins = 40) +
    labs(title = paste("Histogram:", param_name),
         x = param_name, y = "Count") +
    theme_minimal()
  
  # ACF
  acf_obj <- forecast::Acf(df$value, plot = FALSE, lag.max = max_lag)
  acf_df <- tibble(lag = acf_obj$lag[-1], acf = acf_obj$acf[-1])
  p_acf <- ggplot(acf_df, aes(lag, acf)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend = lag, yend = 0)) +
    labs(title = paste("ACF:", param_name),
         x = "Lag", y = "ACF") +
    theme_minimal()
  
  # Geweke
  p_gwk <- geweke_trunc_plot(df$value,
                             param_label = paste(param_name, "district", districtID),
                             max_frac = 0.5, step_frac = 0.05)
  
  list(trace = p_trace, hist = p_hist, acf = p_acf, geweke = p_gwk)
}

diag_plots_Sk <- function(S_arr, districtID, k, max_lag = 50) {
  # S_arr = S_thin
  # districtID = "900959"
  # k=50
  
  str(S_arr)
  df <- extract_Sk(S_arr, districtID, k)
  
  p_trace <- ggplot(df, aes(iter, value)) + geom_line() +
    labs(title = paste("Trace: S[t", k, "] district", districtID),
         x = "Saved draw", y = "S") + theme_minimal()
  
  p_hist <- ggplot(df, aes(value)) + geom_histogram(bins = 40) +
    labs(title = paste("Histogram: S[t", k, "]"),
         x = "S", y = "Count") + theme_minimal()
  
  acf_obj <- forecast::Acf(df$value, plot = FALSE, lag.max = max_lag)
  acf_df <- tibble(lag = acf_obj$lag[-1], acf = acf_obj$acf[-1])
  p_acf <- ggplot(acf_df, aes(lag, acf)) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend = lag, yend = 0)) +
    labs(title = paste("ACF: S[t", k, "]"),
         x = "Lag", y = "ACF") + theme_minimal()
  
  p_gwk <- geweke_trunc_plot(df$value,
                             param_label = paste("Random effect:",k, "district:", districtID),
                             max_frac = 0.5, step_frac = 0.05)
  
  
  list(trace = p_trace, hist = p_hist, acf = p_acf, geweke = p_gwk)
}

plots <- diag_plots_param(P_thin, "102", "logphi")
(plots$trace | plots$hist) / (plots$acf | plots$geweke)
plots <- diag_plots_param(P_thin, "102", "logsigma2")
(plots$trace | plots$hist) / (plots$acf | plots$geweke)
plots <- diag_plots_param(P_thin, "102", "lognu2")
(plots$trace | plots$hist) / (plots$acf | plots$geweke)
plots <- diag_plots_param(P_thin, "102", "sin12")
(plots$trace | plots$hist) / (plots$acf | plots$geweke)

plotsS <- diag_plots_Sk(S_thin, "101", k = 8)
(plotsS$trace | plotsS$hist) / (plotsS$acf | plotsS$geweke)

# --- Posterior fit training set ----
posterior_fit_training <- function(P_draws, S_draws, X_list, y_list,
                                   beta_names = c("b0","sin12","cos12","sin6","cos6"),
                                   district_ids = dimnames(P_draws)[[2]],
                                   probs = c(0.025, 0.5, 0.975)) {
  
  # P_draws= P_thin
  # S_draws= S_thin
  # X_list= X_list
  # y_list= y_list
  # beta_names = c("b0","sin12","cos12","sin6","cos6")
  # district_ids = dimnames(P_draws)[[2]]
  # probs = c(0.025, 0.5, 0.975)
  
  out <- vector("list", length(district_ids))
  post_distr_estimates <- vector("list", length(district_ids))
  names(post_distr_estimates) = district_ids
  names(out) <- district_ids
  
  for (id in district_ids) {
    #id = "101"
    X <- as.matrix(X_list[[id]]) 
    
    y <- as.numeric(y_list[[id]])             # ntime
    ntime <- length(y)
    X = cbind(rep(1,ntime),X)
    ndraw_all <- dim(P_draws)[1]
    p <- ncol(X)
    
    # betas: ndraw x p
    B <- P_draws[, id, beta_names, drop = FALSE]
    B <- matrix(as.numeric(B), nrow = ndraw_all, ncol = p)
    
    # S: ndraw x ntime
    S <- S_draws[, id, , drop = FALSE]
    S <- matrix(as.numeric(S), nrow = ndraw_all, ncol = ntime)
    
    # mu_draws: ndraw x ntime  (fast: B %*% t(X) gives ndraw x ntime)
    mu_draws <- B %*% t(X) + S
    
    logs2 <- as.numeric(P_draws[, id, "logsigma2"])
    lognu <- as.numeric(P_draws[, id, "lognu2"])
    tau2  <- exp(logs2 + lognu)                 # length ndraw
    tau   <- sqrt(tau2)
    
    # posterior predictive draws: yrep = mu + eps, eps ~ N(0, tau2)
    # do it draw-wise efficiently:
    yrep_draws <- mu_draws + matrix(rnorm(ndraw_all * ntime), ndraw_all, ntime) * tau
    
    yrep_mean <- colMeans(yrep_draws)
    yrep_q <- apply(yrep_draws, 2, quantile, probs = probs, names = FALSE)
    rownames(yrep_q) <- c("est.m5.ll","est.m5","est.m5.ul")
    
    estimations = as.data.frame(t(yrep_q))
    out[[id]] = estimations %>% 
      dplyr::mutate(t=t_training, region_id = id)
    
    post_distr_estimates[[id]] = yrep_draws
  }
  
  sum = bind_rows(out)
  return(list(post_distr_estimates,sum))
}

est_results <- posterior_fit_training(P_draws = P_thin, S_draws = S_thin,
                                    X_list = X_list,y_list = y_list,
                                    beta_names = c("b0","sin12","cos12"))
draws_est = est_results[[1]]
est_table = est_results[[2]]
est_table = est_table %>% dplyr::mutate(distcode = as.numeric(region_id))

# Posterior predictive distributions --------------------------------------

posterior_predictive_district <- function(
    districtID, P_draws,S_draws, beta_names = c("b0","sin12","cos12"),
    jitter = 1e-10, seed = 794, t_training, t_test) {
  if (!is.null(seed)) set.seed(seed)
  
  # districtID="101"
  # P_draws = P_thin
  # S_draws = S_thin
  X_test_ind = X_test[[districtID]]
  ntest = dim(X_test_ind)[1]
  X_test_ind = cbind(rep(1,ntest),X_test[[districtID]])
  
  # --- basic checks ---
  stopifnot(districtID %in% dimnames(P_draws)[[2]])
  stopifnot(all(beta_names %in% dimnames(P_draws)[[3]]))
  
  ndraw_all <- dim(P_draws)[1]
  ntime     <- length(t_training)
  ntest     <- length(t_test)
  p         <- ncol(X_test_ind)
  
  draw_idx <- seq_len(ndraw_all)
  
  # pull arrays for this district only (faster inside loop)
  logs2  <- as.numeric(P_draws[draw_idx, districtID, "logsigma2"])
  logphi <- as.numeric(P_draws[draw_idx, districtID, "logphi"])
  lognu2 <- as.numeric(P_draws[draw_idx, districtID, "lognu2"])
  
  B <- P_draws[draw_idx, districtID, beta_names, drop = FALSE]
  B <- matrix(as.numeric(B), nrow = ndraw_all, ncol = p)
  
  S_train <- S_draws[draw_idx, districtID, , drop = FALSE]
  S_train <- matrix(as.numeric(S_train), nrow = ndraw_all, ncol = ntime)
  
  # fixed effects part at test points: ndraw x ntest
  mu_test <- B %*% t(X_test_ind)
  
  # precompute distance matrices once (train/train, train/test, test/test)
  d_tt <- abs(outer(t_training, t_training, "-"))   # ntime x ntime
  d_ts <- abs(outer(t_training, t_test,  "-"))   # ntime x ntest
  d_ss <- abs(outer(t_test,  t_test, "-"))    # ntest x ntest
  
  matern_cor_from_d <- function(d, phi) {
    (1 + d / phi) * exp(-(d / phi))
  }
  
  # output: posterior predictive draws for y_test
  yrep_pred <- matrix(NA_real_, nrow = ndraw_all, ncol = ntest)
  colnames(yrep_pred) <- paste0("pred_t", t_test)
  
  for (k in seq_len(ndraw_all)) {
    #k =1 
    phi  <- exp(logphi[k])
    sig2 <- exp(logs2[k])                      # GP scale
    tau2 <- exp(logs2[k] + lognu2[k])          # obs noise variance
    
    # correlation blocks
    R_tt <- matern_cor_from_d(d_tt, phi)
    R_ts <- matern_cor_from_d(d_ts, phi)
    R_ss <- matern_cor_from_d(d_ss, phi)
    
    # stabilize R_tt
    R_tt <- R_tt + diag(jitter, ntime)
    
    # Cholesky for solves
    cholR <- chol(R_tt)
    
    s_tr <- S_train[k, ]  # length ntime
    
    # conditional mean of s* given s:  C' R^{-1} s
    # where C = R_ts
    Rinvs <- backsolve(cholR, forwardsolve(t(cholR), s_tr))
    m_star <- as.numeric(t(R_ts) %*% Rinvs)  # length ntest
    
    # conditional covariance: sig2 * (R_ss - C' R^{-1} C)
    # compute R^{-1} C efficiently via solves
    Rinvc <- backsolve(cholR, forwardsolve(t(cholR), R_ts))  # ntime x ntest
    V_star <- sig2 * (R_ss - t(R_ts) %*% Rinvc)
    V_star <- V_star + diag(jitter, ntest)
    
    # sample latent s*
    cholV <- chol(V_star)
    s_star <- m_star + as.numeric(t(cholV) %*% rnorm(ntest))
    
    # posterior predictive y* = mu_test + s* + eps*
    eps_star <- rnorm(ntest, mean = 0, sd = sqrt(tau2))
    yrep_pred[k, ] <- mu_test[k, ] + s_star + eps_star
  }
  
  yrep_pred
}

y_pred = vector(mode='list', length = ndistricts)
district_ids <- as.character(distcodes$distcode)

names(y_pred) = district_ids

t_train = training_data %>% ungroup() %>% dplyr::select(t) %>% unique() %>% pull()
t_test = test_data %>% ungroup() %>% dplyr::select(t) %>% unique() %>% pull()

for(id in district_ids){
  y_pred[[id]]<- posterior_predictive_district(
    districtID = id, P_draws = P_thin, S_draws = S_thin, 
    t_training = t_train, t_test = t_test)
}

summarize_ypred <- function(districtID, yrep, probs = c(0.025, 0.5, 0.975)) {
  qs <- apply(yrep, 2, quantile, probs = probs, na.rm = TRUE)
  
  tibble(
    districtID = districtID,
    t = t_test,
    est.m5 = qs[2,],
    est.m5.ll  = qs[1,],
    est.m5.ul = qs[3,])
}

pred_table <- purrr::imap_dfr(
  y_pred,
  ~ summarize_ypred(districtID = .y, yrep = .x)
)
pred_table = pred_table %>% dplyr::rename(distcode = districtID) %>% 
  dplyr::mutate(distcode = as.numeric(distcode))

# Full dataset ------------------------------------------------------------

training_data = training_data %>% left_join(est_table)
test_data = test_data %>% left_join(pred_table)

full_data = ungroup(training_data) %>% 
  add_row(ungroup(test_data))

# --- Plot obs vs est -----
table(full_data$t)
table(full_data$distcode)

full_data %>% 
  dplyr::filter(distcode == 101) %>%
  ggplot(aes(x = t)) + 
  geom_line(aes(y = log_inc, color = "Observed data"), linewidth = 0.5) + 
  geom_line(aes(y = (est.m5), color = "Multivariate GP"), linewidth = 0.5) +
  geom_line(aes(y = (est.m5.ll), color = "Multivariate GP"), linewidth = 0.4, linetype = 2) +
  geom_line(aes(y = (est.m5.ul), color = "Multivariate GP"), linewidth = 0.4, linetype = 2) +
  geom_vline(xintercept = t_training , color = "red", linewidth = 0.5) +
  labs(x = "Time",
       y = "Monthly malaria incidence (log-scale)",
       color = NULL) + 
  scale_color_manual(values = c("Observed data"="#F4538A",
                                "Multivariate GP"="#59D5E0")) +
  theme(
    plot.background = element_rect(fill = rgb(red=41, green=57, blue=69, maxColorValue = 255), color = NA),
    panel.background = element_rect(fill = rgb(red=41, green=57, blue=69, maxColorValue = 255), color = "grey"),
    legend.background = element_rect(fill = rgb(red=41, green=57, blue=69, maxColorValue = 255), color = NA),
    #panel.background = element_rect(fill = rgb(red=41, green=57, blue=69, maxColorValue = 255), color = "grey"),
    plot.title = element_text(color = "white"),
    legend.text = element_text(color = "white"),
    axis.title.x = element_text(color = "white"),
    axis.title.y = element_text(color = "white"),
    axis.text.x = element_text(color = "white", angle=90),
    axis.text.y = element_text(color = "white"),
    axis.line = element_line(linewidth = 0.4, colour="white"), ,
    panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1),  # Very light grey, very thin lines
    panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.1))

# plots  ----- 

full_data_filter = full_data %>% 
  filter(distcode %in% c(101,102,105,114,108,116))

plot_districts = full_data_filter %>% 
  ggplot(aes(x = t)) + 
  geom_line(aes(y = log_inc, color = "Observed data"), linewidth = 0.5) + 
  geom_line(aes(y = (est.m5), color = "Multivariate GP"), linewidth = 0.5) +
  geom_line(aes(y = (est.m5.ll), color = "Multivariate GP"), linewidth = 0.4, linetype = 2) +
  geom_line(aes(y = (est.m5.ul), color = "Multivariate GP"), linewidth = 0.4, linetype = 2) +
  geom_vline(xintercept = t_training , color = "red", linewidth = 0.5) +
  labs(x = "Time",
       y = "Monthly malaria incidence (log-scale)",
       color = NULL) + 
  scale_color_manual(values = c("Observed data"="darkgreen",
                                "Multivariate GP"="#F4538A")) +
  facet_wrap(~distname) + 
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1),
    panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.1)) 

plot_districts
