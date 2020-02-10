// Create a correlation matrix.
//
// get_cor_mat returns a correlation matrix. It is a different parameterization
// for the built-in Stan function cov_exp_quad()
//
// This creates a correlation matrix, \emph{R}, where \deqn{R_{ij} = \prod_{k =
// 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}{R_ij = \prod_{k =
// 1}^{d}\rho_k^{16\left( x_{ik} - x_{jk} \right)^2}}
//
// @param x An n x d matrix, where d is the dimension of the data.
// @param rho A vector of length d. Each value in rho should be between
// 0 and 1.
// @return An n x n correlation matrix

matrix get_cor_mat(matrix x,
                   vector rho) {

  int n = rows(x);
  int d = cols(x);
  real tmp;        // Avoids deep copy warning from Stan
  vector[d] dist4;
  matrix[n, n] R;

  for (i in 1:(n-1)) {
    R[i, i] = 1.0;
    for (j in (i + 1):n) {
      dist4 = 4.0 * (x[i] - x[j])';
      tmp = 1.0;
      for(k in 1:d){
        tmp *= rho[k] ^ (dist4[k]^2);
      }
      R[i, j] = tmp;
      R[j, i] = tmp;
    }
  }
  R[n, n] = 1;
  //return cholesky_decompose(R);
  return R;
}


// Create a covariance matrix.
//
// get_cov_mat_ns returns a covariance matrix (with nugget).
// This creates a covariance matrix, C, where
// \deqn{C =  V^{0.5}RV^{0.5} + sigma^2_\epsilon I}
// where V is a matrix with variances on the diagonal. In the \code{bcgp}
// setting, \emph{R} is a correlation matrix resulting from
// \code{\link{combine_cor_mats}}, \emph{V} is a vector of process variances,
// and \emph{sigma2eps} is the variance of the noise (or nugget).
// @param V A positive vector of length \emph{n}.
// @param R An \emph{n x n} correlation matrix.
// @param sigma2eps A positive scalar representing the variance of the noise
// (or nugget).
// @return An \emph{n x n} covariance matrix

matrix get_cov_mat_ns(vector V,
                      matrix R,
                      real sigma2eps){

  int N = rows(R);
  vector[N] rootV = sqrt(V);
  matrix[N, N] C = quad_form_diag(R, rootV) +
    diag_matrix(rep_vector(sigma2eps,N));

  return C;

}

// Create a covariance matrix.
//
// get_cov_mat_s returns a covariance matrix (with nugget).
// This creates a covariance matrix, C, where
// \deqn{C =  sigma^2 R + sigma^2_\epsilon I}
// In the \code{bcgp}
// setting, \emph{R} is a correlation matrix resulting from
// \code{\link{combine_cor_mats}}, \emph{sigma2} the process variance,
// and \emph{sigma2eps} is the variance of the noise (or nugget).
// @param sigma2 A positive scalar representing the process variance
// @param R An \emph{n x n} correlation matrix.
// @param sigma2eps A positive scalar representing the variance of the noise
// (or nugget).
// @return An \emph{n x n} covariance matrix
matrix get_cov_mat_s(real sigma2,
                     matrix R,
                     real sigma2eps){

  int N = rows(R);
  matrix[N, N] C = sigma2 * R + diag_matrix(rep_vector(sigma2eps,N));

  return C;

}

// Combine correlation matrices.
//
// combine_cor_mats returns a correlation matrix that is a weighted sum of two
// other correlation matrices.
// This creates a correlation matrix, R, where \deqn{R = wG + (1-w)L}. In
// the \code{bcgp} setting, \emph{G} is the global correlation matrix, \emph{L}
// is the local correlation matrix, and \emph{w} is the weight.
//
// @param w A scalar between 0 and 1.
// @param G An n x n correlation matrix, often the result of get_cor_mat
// @param L An n x n correlation matrix, often the result of get_cor_mat
// @return An n x n correlation matrix


matrix combine_cor_mats(real w,
                        matrix G,
                        matrix L){

  int N = rows(G);
  matrix[N, N] R = w*G + (1 - w)*L;
  return R;

}


// Generate a random process conditional on observed.
//
// \code{gp_pred_noncomp_s_rng} returns a vector that is a draw from a GP
// conditional onan already observed vector
//
// @param x Current locations.
// @param y Current observations
// @param x_pred New locations
// @param beta0 Mean of the process
// @param rho Correlation parameter(s) for the process
// @param sigma2eps The random error (or nugget) for the process
// Aparam sigma2 The variance of the process
// @return An n_pred x 1 vector, where n_pred is the number of new locations

vector gp_pred_noncomp_s_rng(matrix x,
                             vector y,
                             matrix x_pred,
                             real beta0,
                             vector rho,
                             real sigma2eps,
                             real sigma2){


  int n_pred = rows(x_pred);
  int n = rows(x);
  int d = cols(x);

  vector[n_pred] y_pred;

  {
    vector[n_pred] y_pred_mean;
    matrix[n_pred, n_pred] y_pred_var_chol;

    int n_all = n_pred + n;
    matrix[n_all, d] x_all = append_row(x_pred, x);
    matrix[n_all, n_all] R_all = get_cor_mat(x_all, rho);
    matrix[n_all, n_all] C_all = get_cov_mat_s(sigma2, R_all, sigma2eps);
    matrix[n_pred, n_pred] C_pred = C_all[1:n_pred, 1:n_pred];
    matrix[n, n] C_L = cholesky_decompose(C_all[(n_pred + 1):, (n_pred + 1):]);
    matrix[n, n_pred] C_cross = C_all[(n_pred + 1):, 1:n_pred];

    vector[n] y_minus_beta0 = y - beta0;
    vector[n] C_L_inv_y_minus_beta0 = mdivide_left_tri_low(C_L, y_minus_beta0);
    matrix[n, n_pred] C_L_inv_C_cross = mdivide_left_tri_low(C_L, C_cross);

    y_pred_mean = beta0 + C_L_inv_C_cross'*C_L_inv_y_minus_beta0;
    y_pred_var_chol =
                  cholesky_decompose(C_pred - C_L_inv_C_cross'*C_L_inv_C_cross);

    y_pred = multi_normal_cholesky_rng(y_pred_mean, y_pred_var_chol);
  }
  return y_pred;
}


// Generate a random process conditional on observed.
//
// \code{gp_pred_comp_ns_rng} returns a vector that is a draw from a GP
// conditional onan already observed vector
//
// @param x Current locations.
// @param y Current observations
// @param x_pred New locations
// @param beta0 Mean of the process
// @param w A scalar between 0 and 1 that is the weight of the global vs. local
// process
// @param rhoL Local correlation parameter(s) for the process
// @param rhoG Global correlation parameter(s) for the process
// @param sigma2eps The random error (or nugget) for the process
// @param V A vector containing the variance of the process at locations in x.
// @param V_pred A vector containing the variance of the process at locations in
// x_pred
// @return An n_pred x 1 vector, where n_pred is the number of new locations
vector gp_pred_comp_ns_rng(matrix x,
                           vector y,
                           matrix x_pred,
                           real beta0,
                           real w,
                           vector rhoL,
                           vector rhoG,
                           real sigma2eps,
                           vector V,
                           vector V_pred){

    int n_pred = rows(x_pred);
    int n = rows(x);
    int d = cols(x);

    vector[n_pred] y_pred;

    {
      vector[n_pred] y_pred_mean;
      matrix[n_pred, n_pred] y_pred_var_chol;

      int n_all = n_pred + n;
      matrix[n_all, d] x_all = append_row(x_pred, x);

      matrix[n_all, n_all] G_all = get_cor_mat(x_all, rhoG);
      matrix[n_all, n_all] L_all = get_cor_mat(x_all, rhoL);
      matrix[n_all, n_all] R_all = combine_cor_mats(w, G_all, L_all);
      vector[n_all] V_all = append_row(V_pred, V);
      matrix[n_all, n_all] C_all = get_cov_mat_ns(V_all, R_all, sigma2eps);

      matrix[n_pred, n_pred] C_pred = C_all[1:n_pred, 1:n_pred];
      matrix[n, n] C_L = cholesky_decompose(C_all[(n_pred + 1):, (n_pred + 1):]);
      matrix[n, n_pred] C_cross = C_all[(n_pred + 1):, 1:n_pred];

      vector[n] y_minus_beta0 = y - beta0;
      vector[n] C_L_inv_y_minus_beta0 = mdivide_left_tri_low(C_L, y_minus_beta0);
      matrix[n, n_pred] C_L_inv_C_cross = mdivide_left_tri_low(C_L, C_cross);

      y_pred_mean = beta0 + C_L_inv_C_cross'*C_L_inv_y_minus_beta0;
      y_pred_var_chol =
                  cholesky_decompose(C_pred - C_L_inv_C_cross'*C_L_inv_C_cross);

      y_pred = multi_normal_cholesky_rng(y_pred_mean, y_pred_var_chol);
    }

    return y_pred;
  }


// gp_pred_noncomp_ns_rng(x, y, x_pred, beta0, rho,
//                                                  sigma2eps, V, V_pred);

// Generate a random process conditional on observed.
//
// \code{gp_pred_noncomp_ns_rng} returns a vector that is a draw from a GP
// conditional onan already observed vector
//
// @param x Current locations.
// @param y Current observations
// @param x_pred New locations
// @param beta0 Mean of the process
// @param rho Correlation parameter(s) for the process
// @param sigma2eps The random error (or nugget) for the process
// @param V A vector containing the variance of the process at locations in x.
// @param V_pred A vector containing the variance of the process at locations in
// x_pred
// @return An n_pred x 1 vector, where n_pred is the number of new locations
vector gp_pred_noncomp_ns_rng(matrix x,
                              vector y,
                              matrix x_pred,
                              real beta0,
                              vector rho,
                              real sigma2eps,
                              vector V,
                              vector V_pred){

    int n_pred = rows(x_pred);
    int n = rows(x);
    int d = cols(x);

    vector[n_pred] y_pred;

    {
      vector[n_pred] y_pred_mean;
      matrix[n_pred, n_pred] y_pred_var_chol;

      int n_all = n_pred + n;
      matrix[n_all, d] x_all = append_row(x_pred, x);

      matrix[n_all, n_all] R_all = get_cor_mat(x_all, rho);
      vector[n_all] V_all = append_row(V_pred, V);
      matrix[n_all, n_all] C_all = get_cov_mat_ns(V_all, R_all, sigma2eps);

      matrix[n_pred, n_pred] C_pred = C_all[1:n_pred, 1:n_pred];
      matrix[n, n] C_L = cholesky_decompose(C_all[(n_pred + 1):, (n_pred + 1):]);
      matrix[n, n_pred] C_cross = C_all[(n_pred + 1):, 1:n_pred];

      vector[n] y_minus_beta0 = y - beta0;
      vector[n] C_L_inv_y_minus_beta0 = mdivide_left_tri_low(C_L, y_minus_beta0);
      matrix[n, n_pred] C_L_inv_C_cross = mdivide_left_tri_low(C_L, C_cross);

      y_pred_mean = beta0 + C_L_inv_C_cross'*C_L_inv_y_minus_beta0;
      y_pred_var_chol =
                  cholesky_decompose(C_pred - C_L_inv_C_cross'*C_L_inv_C_cross);

      y_pred = multi_normal_cholesky_rng(y_pred_mean, y_pred_var_chol);
    }

    return y_pred;
  }


// gp_pred_comp_s_rng(x, y, x_pred, beta0, w, rhoL, rhoG,
//                                              sigma2eps, sigma2);

// Generate a random process conditional on observed.
//
// \code{gp_pred_comp_s_rng} returns a vector that is a draw from a GP
// conditional on an already observed vector
//
// @param x Current locations.
// @param y Current observations
// @param x_pred New locations
// @param beta0 Mean of the process
// @param w A scalar between 0 and 1 that is the weight of the global vs. local
// process
// @param rhoL Local correlation parameter(s) for the process
// @param rhoG Global correlation parameter(s) for the process
// @param sigma2eps The random error (or nugget) for the process
// Aparam sigma2 The variance of the process
// @return An n_pred x 1 vector, where n_pred is the number of new locations

vector gp_pred_comp_s_rng(matrix x,
                          vector y,
                          matrix x_pred,
                          real beta0,
                          real w,
                          vector rhoL,
                          vector rhoG,
                          real sigma2eps,
                          real sigma2){


  int n_pred = rows(x_pred);
  int n = rows(x);
  int d = cols(x);

  vector[n_pred] y_pred;

  {
    vector[n_pred] y_pred_mean;
    matrix[n_pred, n_pred] y_pred_var_chol;

    int n_all = n_pred + n;
    matrix[n_all, d] x_all = append_row(x_pred, x);
    matrix[n_all, n_all] G_all = get_cor_mat(x_all, rhoG);
    matrix[n_all, n_all] L_all = get_cor_mat(x_all, rhoL);
    matrix[n_all, n_all] R_all = combine_cor_mats(w, G_all, L_all);
    matrix[n_all, n_all] C_all = get_cov_mat_s(sigma2, R_all, sigma2eps);
    matrix[n_pred, n_pred] C_pred = C_all[1:n_pred, 1:n_pred];
    matrix[n, n] C_L = cholesky_decompose(C_all[(n_pred + 1):, (n_pred + 1):]);
    matrix[n, n_pred] C_cross = C_all[(n_pred + 1):, 1:n_pred];

    vector[n] y_minus_beta0 = y - beta0;
    vector[n] C_L_inv_y_minus_beta0 = mdivide_left_tri_low(C_L, y_minus_beta0);
    matrix[n, n_pred] C_L_inv_C_cross = mdivide_left_tri_low(C_L, C_cross);

    y_pred_mean = beta0 + C_L_inv_C_cross'*C_L_inv_y_minus_beta0;
    y_pred_var_chol =
                  cholesky_decompose(C_pred - C_L_inv_C_cross'*C_L_inv_C_cross);

    y_pred = multi_normal_cholesky_rng(y_pred_mean, y_pred_var_chol);
  }
  return y_pred;
}
