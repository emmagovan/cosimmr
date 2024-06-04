data {
  int<lower=1> J; // Number of isotopes
  int<lower=1> N; // Number of observations per group
  int<lower=1> K; // Dimension of p and q
  int<lower=1> L1; // Number of covariates in inner
  int<lower=1> L2; //Number of covariates in outer
  matrix[N, J] y; // Data matrix
  matrix[K, J] q; // q matrix
  matrix[K, J] s_mean; // s_mean matrix
  matrix[K, J] c_mean; // c_mean matrix
  matrix[K, J] s_sd; // s_sd matrix
  matrix[K, J] c_sd; // c_sd matrix
  vector[J] sigma_shape; // Prior shape for sigma
  vector[J] sigma_rate; // Prior rate for sigma
  real<lower=0> not_solo; // Adjustment factor for sigma
  matrix[N, L1] X_inner; // Covariates matrix for pack/inside/nested one
    matrix[N, L2] X_outer; // Covariates matrix for pack/inside/nested one
  real<lower=0> hierarchical; // check that data is hierarchical??
}

parameters {
  vector<lower=0>[J] sigma_raw; // log raw sigma values
  matrix[K, L1] beta1; // Matrix of coefficients for inside values (also have intercept here I guess??)
   matrix[K, L2] beta2; // Matrix of coefficients for outside values

}

transformed parameters {
  vector[J] sigma = 0.001 + not_solo * sigma_raw; // Transform sigma_raw
  matrix[N, K] p; // Main parameter
  matrix[N, J] var_y; // Variance for each group
  matrix[N, K] f; // f matrix for CLR prior on p

  for (i in 1:N) {
    for (k in 1:K) {
      f[i,k] = dot_product(X_inner[i,:], beta1[k,:]) + dot_product(X_outer[i,:], beta2[k,:]);
     // f[i, k] = dot_product(X[i,:], beta[k,:]);
     //F is bo + b_region + b_pack (or whatever the nests are)
     // Need to split up X here?? Easier to do it in R i think
    }
  }

  for (i in 1:N) {
    p[i, :] = to_row_vector(softmax(to_vector(f[i, :])));
  }

  for (i in 1:N) {
    for (j in 1:J) {
      var_y[i,j] = dot_product(square(to_vector(p[i,:]) .* q[:, j]), square(s_sd[:, j]) + square(c_sd[:, j]))
                / square(dot_product(to_vector(p[i,:]), q[:, j])) + square(sigma[j]);
    }
  }

}

 model {
  // Prior on betas
  for (k in 1:K) {
    for (l in 1:L1) {
      beta1[k,l] ~ normal(0, 1); // Prior for beta
    }
  }
  
    for (k in 1:K) {
    for (l in 1:L2) {
      beta2[k,l] ~ normal(0, 1); // Prior for beta
    }
  }



  // Prior on sigma_raw
  sigma_raw ~ gamma(sigma_shape, sigma_rate);

  // Likelihood
  for (j in 1:J) {
    for (i in 1:N) {
      real mu_ij = dot_product(to_vector(p[i, :]) .* q[:, j], s_mean[:, j] + c_mean[:, j]) / dot_product(to_vector(p[i, :]), q[:, j]);
      y[i, j] ~ normal(mu_ij, sqrt(var_y[i,j]));
    }
  }
}
