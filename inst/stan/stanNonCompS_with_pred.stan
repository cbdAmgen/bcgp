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
  real<lower = 0> sigma2_alpha;
  real<lower = 0> sigma2_beta;
  real<lower = 0> sigma2eps_alpha;
  real<lower = 0> sigma2eps_beta;

  int<lower=1> n_pred;
  matrix[n_pred, d] x_pred;

}
parameters {

  real beta0;
  vector<lower = 0, upper = 1>[d] rho;
  real<lower=0> sigma2;
  real<lower=0> sigma2eps;

}
model {

  matrix[n, n] C_L;
  vector[n] beta0_vec;
  {

    matrix[n, n] R = get_cor_mat(x, rho);
    matrix[n, n] C = get_cov_mat_s(sigma2, R, sigma2eps);
    C_L = cholesky_decompose(C);
    beta0_vec = rep_vector(beta0, n);

  }

  //beta0 ~ normal(0, 1e6);
  for(i in 1:d)
  rho[i] ~ beta(rho_alpha[i], rho_beta[i]);
  sigma2 ~ gamma(sigma2_alpha, 1/sigma2_beta);
  sigma2eps ~ gamma(sigma2eps_alpha, 1/sigma2eps_beta);

  y ~ multi_normal_cholesky(beta0_vec, C_L);

}
generated quantities{

  vector[n_pred] y_pred = gp_pred_noncomp_s_rng(x, y, x_pred, beta0, rho,
                                                sigma2eps, sigma2);

}
