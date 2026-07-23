d <- abs(outer(inputs$t, inputs$t, "-")) 
cov.mat.GPwosigma = function(phi) {
  matern.cov <- (1 + (d / phi)) * exp(-(d / phi))
  return(matern.cov)
}

# Initialization of the chains ----

coef_names = c("logsigma2","logphi","lognu2","b0",inputs$var_names)

chains_P <- array(NA_real_, dim = c(niter + 1, inputs$ndistricts, length(coef_names)),
                  dimnames = list(NULL, inputs$area_codes, coef_names))

chains_S <- array(NA_real_, dim = c(niter + 1, inputs$ndistricts, inputs$ntimes),
                 dimnames = list(NULL, inputs$area_codes, paste0("t", 1:inputs$ntimes)))

acc.rate = data.frame(matrix(data = NA, nrow = niter-1, ncol = 2*inputs$ndistricts + 1))

# --- init P ---
chains_P[1, , "logsigma2"]  <- inputs$start_values$log.sigma2
chains_P[1, , "logphi"]     <- inputs$start_values$log.phi
chains_P[1, , "lognu2"]     <- inputs$start_values$log.nu2
chains_P[1, , coef_names[4:(inputs$nbetas+3)]] <- 0

# --- init S ---
chains_S[1, , ] <- 0

# Formating Y and X -------------------------------------------------------
str(inputs$training_data)
y_list = lapply(inputs$area_codes, function(id){
  inputs$training_data %>% 
    dplyr::filter(distcode == id) %>% 
    dplyr::select(all_of(inputs$var_y))  %>% dplyr::pull()
})

X_list <- lapply(inputs$area_codes, function(id){
  inputs$X_training %>% 
    dplyr::filter(distcode == id) %>% 
    dplyr::select(-distcode)
})
names(X_list) <- inputs$area_codes
names(y_list) <- inputs$area_codes

# Priors ------------------------------------------------------------------
priorfunction_sigma = function(param, iter_mcmc){
  # param = chains_P
  # iter_mcmc = 1
  param_logsigma2 = param[iter_mcmc,,"logsigma2"]
  
  -0.5 * as.numeric(t(param_logsigma2) %*% inputs$Q %*% param_logsigma2)
}
priorfunction_sigma(chains_P,1)

prior_function_logphi = function(param, iter_mcmc, districtID){
  # param = chains_P
  # iter_mcmc = 1
  # districtID = "322"
  
  distID = as.numeric(districtID)
  
  param_logsigma2 = param[iter_mcmc,districtID,"logsigma2"]
  param_logphi = param[iter_mcmc,districtID,"logphi"]
  
  #Assumed mean for logphi
  mu1 = inputs$start_values %>% 
    dplyr::filter(area_codes == distID) %>% 
    dplyr::select(log.phi) %>% pull(log.phi)
  
  #Assumed mean for logsigma2
  mu2 = inputs$start_values %>% 
    dplyr::filter(area_codes == distID) %>% 
    dplyr::select(log.sigma2) %>% pull(log.sigma2)
  
  #Parts of the matrix
  cov_par = cov_cor_indM_sub[[districtID]]
  cov_sigmaphi = cov_par[2,1]
  var_sigma = cov_par[1,1]
  var_phi = cov_par[2,2]
  
 
  # #Variance of the prior distribution
  #sigma_prior = v_phi
  sigma_prior = max(0.001,var_phi - ((cov_sigmaphi^2)/var_sigma)) #More informative prior
  
  #Mu of the prior distribution
  mu_prior = mu1 + ((cov_sigmaphi)*(1/var_sigma)*(param_logsigma2 - mu2))
  
 dnorm(param_logphi, mean = mu_prior, sd = sqrt(sigma_prior), log = TRUE)
}

prior_function_logphi(chains_P,1,"101")

# Prior distribution for lognu --------------------------------------------

prior_function_lognu = function(param, iter_mcmc, districtID){
  # param = chains_P
  # iter_mcmc = 1
  # districtID = "101"
  distID = as.numeric(districtID)
  
  param_lognu2 = param[iter_mcmc,districtID,"lognu2"]
  
  #Assumed mean for lognu2
  mu_nu = inputs$start_values %>% 
    dplyr::filter(area_codes == distID) %>% 
    dplyr::select(log.nu2) %>% pull(log.nu2)
  
  cov_par = cov_cor_indM_sub[[districtID]]
  se_lognu = sqrt(cov_par[3,3])
  
  dnorm(param_lognu2, mean = mu_nu, sd = se_lognu, log = TRUE)

}

prior_function_lognu(chains_P,1,"101") 


# Prior distribution of the betas -----------------------------------------

prior_function_betas <- function(param, iter_mcmc, districtID) {
  # param = chains_P
  # iter_mcmc = 1
  # districtID = "900959"
  
  param_betas = param[iter_mcmc,districtID,coef_names[4:(3+inputs$nbetas)]]
  S = diag(inputs$nbetas)
  -0.5 * lambda_betas * as.numeric(t(param_betas) %*% S %*% param_betas)
} 

prior_function_betas(chains_P,1,"101")

# Definition of the likelihood  -------------------------------------------

loglik_d <- function(param, param_s, 
                     iter_mcmc, districtID){
  # param = chains_P
  # param_s = chains_S
  # iter_mcmc = 1
  # districtID = "203"
  distID = as.numeric(districtID)
  
  y = y_list[[districtID]]
  
  param_betas = param[iter_mcmc,districtID,coef_names[4:(3+inputs$nbetas)]]
  
  X = as.matrix(cbind(rep(1,inputs$ntimes),
                      X_list[[districtID]])) 
  s = param_s[iter_mcmc, districtID, ]
  
  mu <- (X%*%param_betas) + s
  
  param_lognu2 = param[iter_mcmc,districtID,"lognu2"]
  param_logsigma2 = param[iter_mcmc,districtID,"logsigma2"]
  
  tau2 = exp(param_lognu2+param_logsigma2)
  
  sum(dnorm(y, mean = mu, sd = sqrt(tau2), log = TRUE))
}

loglik_d(chains_P, chains_S,1,"101")

# Prior gaussian process  -------------------------------------------------

gaussian_process<- function(param,param_s,iter_mcmc,districtID){
  # param = chains_P
  # param_s = chains_S
  # iter_mcmc = 1
  # districtID = "900959"
  
  param_logsigma2 = param[iter_mcmc,districtID,"logsigma2"]
  param_logphi = param[iter_mcmc,districtID,"logphi"]
  s = param_s[iter_mcmc, districtID, ]
  
  sigma_gp = cov.mat.GPwosigma(phi = exp(param_logphi)) + diag(1e-10,inputs$ntimes)
  chol_sigma = chol(sigma_gp)
  
  logdet <- inputs$ntimes*param_logsigma2 + 2*sum(log(diag(chol_sigma)))
  
  z <- forwardsolve(t(chol_sigma), s)
  quad <- sum(z^2) * exp(-param_logsigma2)
  
  -0.5*(logdet + quad + inputs$ntimes*log(2*pi))
}

gaussian_process(chains_P,chains_S,1,"101")

# Posterior distributions -------------------------------------------------

post_phi_d <- function(param, param_S, iter_mcmc, districtID){
  prior_function_logphi(param, iter_mcmc, districtID) + 
    gaussian_process(param, param_S, iter_mcmc, districtID)
}
post_phi_d(chains_P, chains_S, 1, "101")

post_for_sigma <- function(param, param_s,iter_mcmc){
  # param = chains_P
  # param_s = chains_S
  # iter_mcmc = 1
  
  lp <- priorfunction_sigma(param, iter_mcmc)  # CAR / joint MVN prior for sigma vector
  dist = as.character(inputs$area_codes)
  
  #dist[19]
  prior = vector(length = length(dist))
  for(i in 1:length(dist)){
    #i = 19
    prior[i] <- prior_function_logphi(param, iter_mcmc, dist[i]) +
      gaussian_process(param, param_s, iter_mcmc, dist[i]) +
      loglik_d(param, param_s, iter_mcmc, dist[i])
  }
  lp+sum(prior)
}
post_for_sigma(chains_P, chains_S, 1)

post_lognu_d <- function(param, param_s, iter_mcmc, districtID){
  prior_function_lognu(param, iter_mcmc, districtID) + loglik_d(param, param_s, iter_mcmc, districtID)
}

# Calculation of the gradients respect to Sigma --------------------------------------------

g_psigma_s = function(param,iter_mcmc){
  param_logsigma2 = as.numeric(param[iter_mcmc,,"logsigma2"])
  g <- as.numeric(- inputs$Q %*% param_logsigma2)
  return(g)
}
g_psigma_s(chains_P,1)

g_pphi_s = function(param, iter_mcmc,districtID){
  # param = chains_P
  # iter_mcmc = 1
  # districtID = "101"
  
  cov_par <- cov_cor_indM_sub[[districtID]]
  logs2  <- param[iter_mcmc, districtID, "logsigma2"]
  logphi <- param[iter_mcmc, districtID, "logphi"]
  
  cov_sigmaphi <- cov_par[2, 1]
  var_sigma    <- cov_par[1, 1]
  var_phi      <- cov_par[2, 2]
  
  a <- cov_sigmaphi / var_sigma
  var_cond <- var_phi - (cov_sigmaphi^2) / var_sigma
  #if (var_cond <= 0) stop("Non-positive conditional variance in phi prior.")
  
  idx <- match(as.numeric(districtID), inputs$start_values$area_codes)
  mu_phi <- inputs$start_values$log.phi[idx]
  mu_sig <- inputs$start_values$log.sigma2[idx]
  mu_cond <- mu_phi + a * (logs2 - mu_sig)
  
  # d/d logsigma2 log N(logphi | mu_cond, var_cond)
  g_phi <- (logphi - mu_cond) * (a / var_cond) 
  
  return(g_phi)
}
g_pphi_s(chains_P,1,"103")

g_pS_s = function(param, param_s, iter_mcmc, districtID, jitter){
  # iter_mcmc = 33
  # param = chains_P
  # param_s = chains_S
  # districtID = 1
  
  logs2  <- param[iter_mcmc, districtID, "logsigma2"]
  logphi <- param[iter_mcmc, districtID, "logphi"]
  s = as.numeric(param_s[iter_mcmc, districtID, ])
  
  sigma_s = cov.mat.GPwosigma(phi = exp(logphi)) + diag(jitter, inputs$ntimes)
  cholR <- chol(sigma_s)
  z <- forwardsolve(t(cholR), s)
  quad <- sum(z^2)  # s' R^{-1} s
  g_gp <- -0.5 * (inputs$ntimes - exp(-logs2) * quad)
  
  return(g_gp)
}
g_pS_s(chains_P,chains_S,1,"101",10^-10)

g_ll_s = function(param, param_s, iter_mcmc, districtID, y, X){
  # iter_mcmc = 1
  # param = chains_P
  # param_s = chains_S
  # districtID = "101"
  
  logs2  <- param[iter_mcmc, districtID, "logsigma2"]
  logphi <- param[iter_mcmc, districtID, "logphi"]
  s = param_s[iter_mcmc, districtID, ]
  param_betas = param[iter_mcmc,districtID,4:(3+inputs$nbetas)]
  param_lognu2 = param[iter_mcmc,districtID,"lognu2"]
  param_logsigma2 = param[iter_mcmc,districtID,"logsigma2"]
  
  tau2 = exp(param_lognu2+param_logsigma2)
  
  mu <- (X%*%param_betas) + s
  
  rss  <- sum((y - mu)^2)
  g_lik <- 0.5 * (rss / tau2 - inputs$ntimes)
  return(g_lik)
}

total_gradient_sigma = function(param, param_S, iter_mcmc){
  gradient_sigma = as.vector(0)
  for(i in 1:inputs$ndistricts){
    districtID = as.character(inputs$area_codes)[i]
    y = y_list[[districtID]]
    
    X = as.matrix(cbind(rep(1,inputs$ntimes),X_list[[districtID]])) 
    
    gradient_sigma[i] = g_pphi_s(param, iter_mcmc, districtID) + 
      g_pS_s(param, param_S, iter_mcmc, districtID, jitter = 10^(-10)) + 
      g_ll_s(param, param_S, iter_mcmc, districtID, y, X)
  }
  gradient_sigma = gradient_sigma + g_psigma_s(param, iter_mcmc)
  return(gradient_sigma)
}
total_gradient_sigma(chains_P, chains_S,1)


# Calculation of the gradients respect to phi -----------------------------

g_pphi_p = function(param, iter_mcmc, districtID){
  # district row
  idx <- match(as.numeric(districtID), inputs$area_codes)
  if (is.na(idx)) stop("districtID not found in cov_pars")
  
  logs2  <- as.numeric(param[iter_mcmc, districtID, "logsigma2"])
  logphi <- as.numeric(param[iter_mcmc, districtID, "logphi"])
  
  cov_par = cov_cor_indM_sub[[districtID]]
  cov_sigmaphi <- cov_par[2,1]
  var_sigma    <- cov_par[1,1]
  var_phi      <- cov_par[2,2]
  
  a <- cov_sigmaphi / var_sigma
  var_cond <- var_phi - (cov_sigmaphi^2)/var_sigma
  
  if (var_cond <= 0) stop("Non-positive conditional variance in phi prior.")
  
  mu_phi <- cov.pars$logphi[idx]
  mu_sig <- cov.pars$logsigma[idx]   # adjust if your column is named differently
  mu_cond <- mu_phi + a*(logs2 - mu_sig)
  
  # d/d logphi log N(logphi | mu_cond, var_cond)
  gradient = - (logphi - mu_cond) / var_cond
  return(gradient)
}
g_pphi_p(chains_P,1,"101")

g_pS_p = function(param, param_s, iter_mcmc, districtID, jitter){
  logs2  <- as.numeric(param[iter_mcmc, districtID, "logsigma2"])
  logphi <- as.numeric(param[iter_mcmc, districtID, "logphi"])
  s   <- as.numeric(param_s[iter_mcmc, districtID, ])
  
  phi <- exp(logphi)
  R <- cov.mat.GPwosigma(phi) + diag(10^(-10), inputs$ntimes)
  cholR <- chol(R)
  
  # dR/dphi for R=(1+d/phi)exp(-d/phi)
  # uses global d (distance matrix)
  dR_dphi <- (d^2 / phi^3) * exp(-(d / phi))
  
  # R^{-1} s
  Rinvs <- backsolve(cholR, forwardsolve(t(cholR), s))
  
  # R^{-1} dR
  RinvdR <- backsolve(cholR, forwardsolve(t(cholR), dR_dphi))
  
  trace_term <- sum(diag(RinvdR))
  quad_term  <- as.numeric(t(Rinvs) %*% dR_dphi %*% Rinvs)
  
  dlog_dphi <- -0.5 * trace_term + 0.5 * exp(-logs2) * quad_term
  
  # chain rule to logphi
  gradient = phi * dlog_dphi
  return(gradient)
}
g_pS_p(chains_P,chains_S, 1,"101",10^(-10))

grad_logphi_total <- function(param, param_S, iter_mcmc, districtID, jitter = 1e-10){
  gradient = g_pphi_p(param, iter_mcmc, districtID) +
    g_pS_p(param, param_S, iter_mcmc, districtID, jitter)
  return(gradient)
}

grad_logphi_total(chains_P,chains_S, 1,"101")

# Proposal value log(sigma2) ----------------------------------------------

proposalfunctionlogsigmaMulti <- function(param, param_S, iter_mcmc, eps){
  x <- as.numeric(param[iter_mcmc, , "logsigma2"])     
  g <- total_gradient_sigma(param, param_S, iter_mcmc)                          # numeric length D
  
  mu_prop <- x + 0.5 * eps^2 * g
  x_new <- as.numeric(mvtnorm::rmvnorm(1, mean = mu_prop, sigma = eps^2 * diag(length(x))))
  x_new
}
proposalfunctionlogsigmaMulti(chains_P,chains_S,1,0.01)

# Proposal value log(phi) and log(nu2) ------------------------------------

proposal_logphi_mala <- function(param, param_S, iter_mcmc, districtID, eps, jitter = 1e-10){
  # Current value
  x_curr <- as.numeric(param[iter_mcmc, districtID, "logphi"])
  
  # Gradient of log posterior wrt logphi at current value
  g_curr <- as.numeric(grad_logphi_total(param, param_S, iter_mcmc, districtID, jitter))
  
  # MALA drift mean
  mu_prop <- x_curr + 0.5 * eps^2 * g_curr
  
  # Sample proposal: N(mu_prop, eps^2)
  x_prop <- rnorm(1, mean = mu_prop, sd = eps)
  
  return(x_prop)
}

proposal_lognu_rw <- function(param, iter_mcmc, districtID){
  idx <- match(as.numeric(districtID), inputs$area_codes)
  x <- as.numeric(param[iter_mcmc, districtID, "lognu2"])
  sd_lognu = sqrt(cov_cor_indM_sub[[districtID]][3,3])
  nu_prop = rnorm(1, mean = x, sd = sd_lognu)
  return(nu_prop)
}

# Function for epsilon ----------------------------------------------------

hfun <- function(x) {
  epsil <- 1e-04
  Apar <- 10 ^ 4
  if (x < epsil) {val <- epsil} 
  else if (x >= epsil && x <= Apar) {val <- x} 
  else{val <- Apar}
  return(val)
}

# Quotient for ar step ----------------------------------------------------

log_qratio_logsigma <- function(param, param_S, iter_mcmc, eps){
  logsigma_curr <- as.numeric(param[iter_mcmc,     , "logsigma2"])
  logsigma_prop <- as.numeric(param[iter_mcmc + 1, , "logsigma2"])
  
  # forward mean mu_fwd = x_curr + 0.5 eps^2 g(curr)
  g_curr <- as.numeric(total_gradient_sigma(param, param_S, iter_mcmc))
  mu_fwd <- logsigma_curr + 0.5 * eps^2 * g_curr
  
  # backward mean mu_bwd = x_prop + 0.5 eps^2 g(prop)
  g_prop <- as.numeric(total_gradient_sigma(param, param_S, iter_mcmc+1))
  mu_bwd <- logsigma_prop + 0.5 * eps^2 * g_prop
  
  # covariance = eps^2 I  => constants cancel, but easiest is mvtnorm
  Sigma <- eps^2 * diag(inputs$ndistricts)
  
  log_q_bwd <- mvtnorm::dmvnorm(logsigma_curr, mean = mu_bwd, sigma = Sigma, log = TRUE)
  log_q_fwd <- mvtnorm::dmvnorm(logsigma_prop, mean = mu_fwd, sigma = Sigma, log = TRUE)
  
  log_q_bwd - log_q_fwd
}

log_qratio_logphi <- function(param, param_S, iter_mcmc, districtID, eps, jitter = 1e-10, grad_logphi_total){
  logphi_curr <- as.numeric(param[iter_mcmc,     districtID, "logphi"])
  logphi_prop <- as.numeric(param[iter_mcmc + 1, districtID, "logphi"])  # assumes proposal stored
  
  g_curr <- as.numeric(grad_logphi_total(param, param_S, iter_mcmc, districtID, jitter))
  mu_fwd <- logphi_curr + 0.5 * eps^2 * g_curr
  
  g_prop <- as.numeric(grad_logphi_total(param, param_S, iter_mcmc + 1, districtID, jitter))
  mu_bwd <- logphi_prop + 0.5 * eps^2 * g_prop
  
  log_q_bwd <- dnorm(logphi_curr, mean = mu_bwd, sd = eps, log = TRUE)
  log_q_fwd <- dnorm(logphi_prop, mean = mu_fwd, sd = eps, log = TRUE)
  
  log_q_bwd - log_q_fwd
}

# Step-size storage -------------------------------------------------------

eps_sigma <- numeric(niter)
eps_sigma[1] <- epsilon_0

eps_phi <- matrix(epsilon_0, nrow = niter, ncol = inputs$ndistricts,
                  dimnames = list(NULL, as.character(inputs$area_codes)))

# Function updates district-specific --------------------------------------

update_one_district <- function(param, param_s, iter_mcmc, districtID, eps_phi_j){
  # param = chains_P
  # param_s = chains_S
  # iter_mcmc = 1
  # districtID = "103"
  # eps_phi_j = epsilon_0
  # 
  # start from current state at j
  param_all_next <- param[iter_mcmc, districtID, ]
  random_next <- param_s[iter_mcmc, districtID, ]
  
  ## Step 1: Sampling log_phi
  # Proposal value
  prop_phi <- proposal_logphi_mala(param, param_s, iter_mcmc, districtID, 
                                   eps = eps_phi_j)
  
  # Store the proposal into a temp "j+1" slot for gradient/qratio eval
  Ptmp <- param
  Ptmp[iter_mcmc+1, districtID, ] <- param[iter_mcmc, districtID, ]
  Ptmp[iter_mcmc+1, districtID, "logphi"] <- prop_phi
  Sarr_tmp <- param_s
  Sarr_tmp[iter_mcmc+1, districtID, ] <- param_s[iter_mcmc, districtID, ]
  
  # log posterior for phi block (phi prior + GP prior for S)
  lp_curr <- post_phi_d(param, param_s, iter_mcmc, districtID)
  lp_prop <- post_phi_d(Ptmp, Sarr_tmp, iter_mcmc+1, districtID)  
  
  # q-ratio
  log_q <- log_qratio_logphi(Ptmp, Sarr_tmp, iter_mcmc, districtID, 
                             eps = eps_phi_j, jitter = 1e-10, 
                             grad_logphi_total = grad_logphi_total)
  
  # log of the acceptance prob
  log_alpha <- (lp_prop - lp_curr) + log_q
  # AR step 
  acc_phi <- as.integer(log(runif(1)) < log_alpha)
  
  if (acc_phi == 1L) param_all_next["logphi"] <- prop_phi
  
  ## Step 2: Sampling lognu2 
  # Proposal value
  prop_lognu = proposal_lognu_rw(param, iter_mcmc, districtID)
  
  # Store the proposal into a temp "j+1" slot for gradient/qratio eval
  Ptmp <- param
  Ptmp[iter_mcmc+1, districtID, ] <- param[iter_mcmc, districtID, ]
  Ptmp[iter_mcmc+1, districtID, "lognu2"] <- prop_lognu
  Sarr_tmp <- param_s
  Sarr_tmp[iter_mcmc+1, districtID, ] <- param_s[iter_mcmc, districtID, ]
  
  # log posterior for lognu 
  lp_curr_nu = post_lognu_d(param, param_s,iter_mcmc, districtID)
  lp_prop_nu = post_lognu_d(Ptmp, Sarr_tmp,iter_mcmc+1, districtID)
  
  # Acceptance rejection step 
  log_alpha_nu = lp_prop_nu - lp_curr_nu
  acc_nu = as.integer(log(runif(1)) < log_alpha_nu)
  if (acc_nu == 1L) param_all_next["lognu2"] <- prop_lognu
  
  ## Step 3: Sampling betas given sigma, nu S, using the penalty prior on the betas
  # Use your existing formula if you like; but avoid solve(solve()) if possible.
  Xd <- as.matrix(cbind(rep(1,inputs$ntimes),X_list[[districtID]]))
  yd <- y_list[[districtID]]
  tau2 <- exp(param_all_next["lognu2"] + param_all_next["logsigma2"])
  lam <- lambda_betas
  
  # posterior precision: (1/tau2) X'X + lam S
  S = diag(inputs$nbetas)
  Prec <- ((crossprod(Xd) / tau2) +  (lam*S)) + diag(10^-(10),inputs$nbetas)
  cholP <- chol(Prec)
  
  # mean: Prec^{-1} (1/tau2) X'(y - s)
  rhs <- crossprod(Xd, yd - random_next) / tau2
  mu  <- backsolve(cholP, forwardsolve(t(cholP), rhs))
  
  # sample: mu + Prec^{-1/2} z
  z <- rnorm(ncol(Xd))
  sample_b <- as.numeric(mu + backsolve(cholP, z))
  param_all_next[4:(3+inputs$nbetas)] <- sample_b
  
  ## Step 4: Sampling the random effects 
  phi <- exp(param_all_next["logphi"])
  sigma = exp(param_all_next["logsigma2"]) 
  
  R <- cov.mat.GPwosigma(phi) + diag(1e-10, inputs$ntimes)
  
  # Posterior precision for S: (1/tau2) I + exp(-logsigma2) R^{-1}
  cholR <- chol(R)
  # compute R^{-1} as operator via solves
  # build A = (1/tau2)I + exp(-logs2)R^{-1}
  # We'll apply via dense matrix because ntime might be ~100; ok.
  Rinvd <- backsolve(cholR, forwardsolve(t(cholR), diag(inputs$ntimes))) # could be heavy; if ntime small OK
  A <- diag(1/tau2, inputs$ntimes) + (1/sigma)* Rinvd
  cholA <- chol(A)
  
  muS_rhs <- (yd - as.numeric(Xd %*% sample_b)) / tau2
  muS <- backsolve(cholA, forwardsolve(t(cholA), muS_rhs))
  z_re = rnorm(inputs$ntimes)
  random_next <- as.numeric(muS + backsolve(cholA, z_re))
  
  list(id = districtID,Pnext = param_all_next,Snext = random_next,
       acc_phi = acc_phi,acc_nu = acc_nu)
}

# Sampler -----------------------------------------------------------------

sampler_mcmc_array <- function(niter){
  #acc_s = vector(0)
  # niter = 100
  ids <- as.character(inputs$area_codes)
  
  for (j in 1:(niter-1)) {
    #j =1
    
    ## Step 1: Sampling logsigma 
    chains_P[j+1,,] <- chains_P[j,,]
    chains_S[j+1,,]  <- chains_S[j,,]
    
    prop_sigma <- proposalfunctionlogsigmaMulti(chains_P, chains_S, j, eps_sigma[j])
    chains_P[j+1, , "logsigma2"] <- prop_sigma
    
    # MH for sigma block
    lp_curr <- post_for_sigma(chains_P, chains_S, j)
    lp_prop <- post_for_sigma(chains_P, chains_S, j+1)
    
    log_q <- log_qratio_logsigma(chains_P, chains_S, j, eps_sigma[j])
    
    log_alpha <- (lp_prop - lp_curr) + log_q
    acc_s <- as.integer(log(runif(1)) < log_alpha)
  
    if (acc_s == 0L) {
      chains_P[j+1, , "logsigma2"] <- chains_P[j, , "logsigma2"]
    }
    
    # adapt eps_sigma 
    acc_prob <- min(1, exp(log_alpha))
    eps_sigma[j+1] <- (hfun(sqrt(eps_sigma[j]) + (1/(j+1))*(acc_prob - 0.57)))^2
    
    ## --- (2) district updates in parallel ---
    # make sure chains_P[j+1, , ] has updated sigma already
    # and everything else carried forward.
    
    # PARALLEL (Unix): parallel::mclapply
    # Windows: use parallel::parLapply with cluster or future.apply
    res_list <- parallel::mclapply(
      ids,
      function(id){
        update_one_district(chains_P, chains_S, j+1, id, eps_phi[j,id])},
      #mc.cores = max(1L, parallel::detectCores() - 1L)
      mc.cores = max(1L, 7 - 1L)
    )
    
    # write back results
    for (k in seq_along(res_list)) {
      #k=1
      r <- res_list[[k]]
      id <- r$id
      chains_P[j+1, id, ] <- r$Pnext
      chains_S[j+1, id, ]  <- r$Snext
      
      # adapt eps_phi for this district
      eps_phi[j+1, id] <- (hfun(sqrt(eps_phi[j, id]) + (1/(j+1))*(r$acc_phi - 0.57)))^2
      # track acc if you want
    }
    print(paste("iter",j,"done"))
  }
  list(P = chains_P, S = chains_S, eps_sigma = eps_sigma, eps_phi = eps_phi)
}


