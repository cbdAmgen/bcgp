functions {

#include /functions/common_functions.stan

}
data {

  int<lower=1> n;
  int<lower=1> d;
  vector[n] y;
  matrix[n, d] x;

  vector<lower = 0>[d] rho_alpha;
  vector<lower = 0>[d] rho_beta;
  real muV_betaV;
  real<lower = 0> muV_sigma2;
  vector<lower = 0>[d] rhoV_alpha;
  vector<lower = 0>[d] rhoV_beta;
  real<lower = 0> sigma2V_alpha;
  real<lower = 0> sigma2V_beta;
  real<lower = 0> sigma2eps_alpha;
  real<lower = 0> sigma2eps_beta;

}
parameters {

  real beta0;
  vector<lower = 0, upper = 1>[d] rho;
  real muV;
  real<lower = 0> sigma2V;
  vector<lower = 0, upper = 1>[d] rhoV;
  vector[n] logV_raw;
  real<lower=0> sigma2eps;

}
transformed parameters{

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

    matrix[n, n] R = get_cor_mat(x, rho);
    matrix[n, n] C = get_cov_mat_ns(V, R, sigma2eps);
    C_L = cholesky_decompose(C);
    beta0_vec = rep_vector(beta0, n);

  }

  for(i in 1:d){
    rho[i] ~ beta(rho_alpha[i], rho_beta[i]);
    rhoV[i] ~ beta(rhoV_alpha[i], rhoV_beta[i]);
  }
  sigma2eps ~ gamma(sigma2eps_alpha, 1/sigma2eps_beta);

  muV ~ normal(muV_betaV, muV_sigma2);
  sigma2V ~ inv_gamma(sigma2V_alpha, 1/sigma2V_beta);
  logV_raw ~ normal(0, 1);

  y ~ multi_normal_cholesky(beta0_vec, C_L);

}
