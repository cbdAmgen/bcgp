functions {

#include /functions/common_functions.stan

}
data {

  int<lower=1> n;
  int<lower=1> d;
  vector[n] y;
  matrix[n, d] x;

  real<lower = 0, upper = 1> w_lower;
  real<lower = w_lower, upper = 1> w_upper;
  real<lower = 0> w_alpha;
  real<lower = 0> w_beta;
  vector<lower = 0>[d] rhoG_alpha;
  vector<lower = 0>[d] rhoG_beta;
  vector<lower = 0>[d] rhoL_alpha;
  vector<lower = 0>[d] rhoL_beta;
  real muV_betaV;
  real<lower = 0> muV_sigma2;
  vector<lower = 0>[d] rhoV_alpha;
  vector<lower = 0>[d] rhoV_beta;
  real<lower = 0> sigma2V_alpha;
  real<lower = 0> sigma2V_beta;
  real<lower = 0> sigma2eps_alpha;
  real<lower = 0> sigma2eps_beta;

  int<lower=1> n_pred;
  matrix[n_pred, d] x_pred;

}
parameters {

  real beta0;
  real<lower = 0, upper = 1> w_raw;
  vector<lower = 0, upper = 1>[d] rhoG;
  vector<lower = 0, upper = 1>[d] rhoL_raw;
  real muV;
  real<lower = 0> sigma2V;
  vector<lower = 0, upper = 1>[d] rhoV;
  vector[n] logV_raw;
  real<lower=0> sigma2eps;

}
transformed parameters{

  real<lower = w_lower, upper = w_upper> w = w_lower +
                                             (w_upper - w_lower)*w_raw;
  vector<lower = 0, upper = 1>[d] rhoL = rhoG .* rhoL_raw;
  vector[n] V;

  {

    matrix[n, n] R_K = get_cor_mat(x, rhoV);
    matrix[n, n] K = get_cov_mat_s(sigma2V, R_K, 1e-10);
    matrix[n, n] K_L = cholesky_decompose(K);
    vector[n] muV_vec = rep_vector(muV, n);
    vector[n] logV = muV_vec + K_L*logV_raw;
    V = exp(logV);

  }

}
model {

  matrix[n, n] C_L;
  vector[n] beta0_vec;

  {

    matrix[n, n] G = get_cor_mat(x, rhoG);
    matrix[n, n] L = get_cor_mat(x, rhoL);
    matrix[n, n] R = combine_cor_mats(w, G, L);
    matrix[n, n] C = get_cov_mat_ns(V, R, sigma2eps);
    C_L = cholesky_decompose(C);
    beta0_vec = rep_vector(beta0, n);

  }


  w_raw ~ beta(w_alpha, w_beta);
  for(i in 1:d){
    rhoG[i] ~ beta(rhoG_alpha[i], rhoG_beta[i]);
    rhoL_raw[i] ~ beta(rhoL_alpha[i], rhoL_beta[i]);
    rhoV[i] ~ beta(rhoV_alpha[i], rhoV_beta[i]);
  }
  sigma2eps ~ gamma(sigma2eps_alpha, 1/sigma2eps_beta);

  muV ~ normal(muV_betaV, muV_sigma2);
  sigma2V ~ inv_gamma(sigma2V_alpha, 1/sigma2V_beta);
  logV_raw ~ normal(0, 1);

  y ~ multi_normal_cholesky(beta0_vec, C_L);

}
generated quantities{

  vector[n_pred] V_pred = exp(gp_pred_noncomp_s_rng(x, log(V), x_pred, muV, rhoV,
                                                1e-8, sigma2V));
  vector[n_pred] y_pred = gp_pred_comp_ns_rng(x, y, x_pred, beta0, w, rhoL,
                                              rhoG, sigma2eps, V, V_pred);

}

