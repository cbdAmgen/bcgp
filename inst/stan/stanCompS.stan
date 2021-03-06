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
  real<lower = 0> sigma2_alpha;
  real<lower = 0> sigma2_beta;
  real<lower = 0> sigma2eps_alpha;
  real<lower = 0> sigma2eps_beta;

}
parameters {

  real beta0;
  real<lower = 0, upper = 1> w_raw;
  vector<lower = 0, upper = 1>[d] rhoG;
  vector<lower = 0, upper = 1>[d] rhoL_raw;
  real<lower=0> sigma2;
  real<lower=0> sigma2eps;

}
transformed parameters{

  real<lower = w_lower, upper = w_upper> w = w_lower +
                                             (w_upper - w_lower)*w_raw;
  vector<lower = 0, upper = 1>[d] rhoL = rhoG .* rhoL_raw;

}
model {

  matrix[n, n] C_L;
  vector[n] beta0_vec;
  {

    matrix[n, n] G = get_cor_mat(x, rhoG);
    matrix[n, n] L = get_cor_mat(x, rhoL);
    matrix[n, n] R = combine_cor_mats(w, G, L);
    matrix[n, n] C = get_cov_mat_s(sigma2, R, sigma2eps);
    C_L = cholesky_decompose(C);
    beta0_vec = rep_vector(beta0, n);

  }

  w_raw ~ beta(w_alpha, w_beta);
  for(i in 1:d){
    rhoG[i] ~ beta(rhoG_alpha[i], rhoG_beta[i]);
    rhoL_raw[i] ~ beta(rhoL_alpha[i], rhoL_beta[i]);
  }
  sigma2 ~ gamma(sigma2_alpha, 1/sigma2_beta);
  sigma2eps ~ gamma(sigma2eps_alpha, 1/sigma2eps_beta);

  y ~ multi_normal_cholesky(beta0_vec, C_L);

}
