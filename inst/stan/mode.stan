functions {
   // assumes A is SPD
    matrix square_root(matrix A, int R) {
        matrix[R, R] eigvecs = eigenvectors_sym(A);
        vector[R] eigvals = eigenvalues_sym(A);
        return eigvecs * diag_matrix(sqrt(eigvals)) / eigvecs;
    }

    real log_norm_const(int d, real rho) {
      real num = (d-1)*log(2) + log(1+inv_sqrt(1-pow(rho, 2)));
      real den = d * log(sqrt(1 - rho) + sqrt(1 + rho));
      return num - den;
    }

    real mymed_lpdf(vector x, vector m, matrix sigma_inv, matrix sigma_inv_sqrt, vector nu, real r, int d) {
      real out = multi_normal_prec_lpdf(x | m, sigma_inv);
      out -= (log_norm_const(d, r));
      out -= r / 2* (x-m)'* sigma_inv_sqrt * nu * sqrt(quad_form(sigma_inv, x-m));
      return out;
    }
}

data {
  int<lower=1> d;         // data dimension
  int<lower=0> n;         // number of observations
  vector[d] Y[n];         // observations
}

parameters {
  vector[d] mu;                     // mean
  real<lower=0,upper=1> rho;        // asymmetry parameter
  cov_matrix[d] Sigma_inv;          // precision matrix
  vector<lower=0>[d] diag_sigma_sq; // mu hyperprior cov mat
  unit_vector[d] nu_tilde;          // direction (tilde version)
}

transformed parameters {
  matrix[d, d] sigma_inv_sqrt = square_root(Sigma_inv, d);
}


model {
  mu ~ normal(rep_vector(0, d), diag_sigma_sq);
  for (j in 1:d) {
    diag_sigma_sq[j] ~ inv_gamma(1, 1);
  }
  rho ~ uniform(0, 1);
  Sigma_inv ~ wishart(d, diag_matrix(rep_vector(1, d)));
  for (i in 1:n) {
    Y[i] ~ mymed( mu, Sigma_inv, sigma_inv_sqrt, nu_tilde, rho, d);
  }
}
