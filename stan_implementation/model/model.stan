functions {
  vector GP_1D(array[] real x, real rho, real alpha, vector eta) {
    return cholesky_decompose(
      add_diag(gp_matern32_cov(x, alpha, rho), 1e-9)
    ) * eta;
  }

  real partial_sum(
    array[] int site_ids,
    int start, int end,
    array[] vector y,
    array[] matrix seasonality,
    array[] real time,
    array[] real intercept,
    array[] real log_phi,
    vector log_sigma,
    array[] vector eta,
    array[] vector beta_sincos,
    array[] real log_tau,
    real mu_log_phi,
    real gamma_phi_sigma,
    real sd_log_phi
  ) {
    real lp = 0;

    for (n in 1:size(site_ids)) {
      int i = site_ids[n];

      lp += normal_lpdf(intercept[i] | 0, 3);

      lp += normal_lpdf(
        log_phi[i] |
        mu_log_phi + gamma_phi_sigma * log_sigma[i],
        sd_log_phi
      );

      lp += normal_lpdf(log_tau[i] | log(0.5), 1);
      lp += std_normal_lpdf(eta[i]);
      lp += normal_lpdf(beta_sincos[i] | 0, 3);

      lp += normal_lpdf(
        y[i] |
        intercept[i] +
        seasonality[i] * beta_sincos[i] +
        GP_1D(time, exp(log_phi[i]), exp(log_sigma[i]), eta[i]),
        exp(log_tau[i])
      );
    }

    return lp;
  }
}

data {
  int<lower=1> n_times;
  int<lower=1> n_sites;
  int<lower=1> n_sincos;
  int<lower=0> n_edges;

  array[n_sites] vector[n_times] y;
  array[n_times] real time;
  array[n_sites] matrix[n_times, n_sincos] seasonality;

  array[n_edges] int<lower=1, upper=n_sites> node1;
  array[n_edges] int<lower=1, upper=n_sites> node2;

  int<lower=1> grainsize;
}

transformed data {
  array[n_sites] int site_ids;

  for (i in 1:n_sites) {
    site_ids[i] = i;
  }
}

parameters {
  array[n_sites] real intercept;
  array[n_sites] real log_phi;
  vector[n_sites] log_sigma;
  array[n_sites] vector[n_times] eta;
  array[n_sites] vector[n_sincos] beta_sincos;
  array[n_sites] real log_tau;

  real<lower=0> sigma_icar;

  real mu_log_phi;
  real gamma_phi_sigma;
  real<lower=0> sd_log_phi;
}

model {
  real prec_sigmas = inv_square(sigma_icar);

  sigma_icar ~ normal(0, 0.3);

  target += -0.5 * prec_sigmas *
    dot_self(log_sigma[node1] - log_sigma[node2])
    - (n_sites - 1) * log(sigma_icar);

  sum(log_sigma) ~ normal(0, 0.01 * n_sites);

  mu_log_phi ~ normal(log(0.35), 0.4);
  gamma_phi_sigma ~ normal(0, 0.5);
  sd_log_phi ~ normal(0.5, 0.2);

  target += reduce_sum(
    partial_sum,
    site_ids,
    grainsize,
    y,
    seasonality,
    time,
    intercept,
    log_phi,
    log_sigma,
    eta,
    beta_sincos,
    log_tau,
    mu_log_phi,
    gamma_phi_sigma,
    sd_log_phi
  );
}

generated quantities {
  array[n_sites] vector[n_times] y_rep;
  array[n_sites] vector[n_times] mu;
  array[n_sites, n_times] real log_lik;

  for (i in 1:n_sites) {
    vector[n_times] temporal_effect;

    temporal_effect = GP_1D(
      time,
      exp(log_phi[i]),
      exp(log_sigma[i]),
      eta[i]
    );

    mu[i] =
      intercept[i] +
      seasonality[i] * beta_sincos[i] +
      temporal_effect;

    for (t in 1:n_times) {
      y_rep[i][t] = normal_rng(mu[i][t], exp(log_tau[i]));
      log_lik[i, t] = normal_lpdf(y[i][t] | mu[i][t], exp(log_tau[i]));
    }
  }
}
