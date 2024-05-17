#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]




// [[Rcpp::export]]
arma::mat rMVNormCpp_pxr(int n, arma::vec Mean, arma::mat Var) {
  int ncols = Var.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(Mean, 1, n).t() + Y * arma::chol(Var);
}

//The part that is called epsilon in the algorithm write up - I'm renamming
// kappa because we have several epsilons already and its getting confusing

//And now the part that is called epsilon in mixsiar maths is called omicron
//why does everyone like epsilon so much
//So now we have y ~ N(sum(p*q*mu), sum(p*q*sigma_TDF) x omicron + sigma^2)


//simulates S samples of theta - so several rnorms + rnorm for log(sigma^2) + rnorm for log(omicron^2)
//[[Rcpp::export]]
arma::mat sim_thetacpp_pxr(int S, arma::vec lambda, int n_sources,
                       int n_tracers, int n_cov, bool solo,
                       arma::mat kappa){
  
  int ncnsnt = n_sources * n_cov + (n_tracers * 2);
  
  arma::mat theta(S, (ncnsnt));
  
  arma::vec mean = lambda.subvec(0, ncnsnt - 1);
  
  int mat_size = ((ncnsnt) * (ncnsnt+1))/2;
  
  arma::vec sig = lambda.subvec(ncnsnt, ncnsnt + mat_size -1);
  
  
  arma::mat chol_var(ncnsnt, ncnsnt, arma::fill::zeros);
  
  int count = 0;
  
  for(int i = 0; i<ncnsnt; i++){
    for(int j = 0; j<ncnsnt; j++){
      if (j <= i){
        count +=1;
        chol_var((i),(j)) = sig(count-1);
        
      }
      
    }
  }

  
  arma::mat Leps = chol_var*kappa;
  
  for(int i=0; i<S; i++){
    theta.row(i) = (mean + Leps.col(i)).t();
  }
  
  
  return theta;
}



// calculates p - which is exp(f) / sum exp(f)
// and f is X * B
// X is a matrix and B is a vector

//[[Rcpp::export]]
arma::mat hfn_pxr(arma::vec theta,
              int n_sources, int n, int n_cov,
              arma::mat x_scaled){
  
  arma::mat beta(n_cov, n_sources);
  
  for(int i = 0; i<n_cov; i++){
    for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
  }
  
  arma::mat f = x_scaled * beta;
  arma::mat expf = exp(f);
  
  
  arma::vec sumexpf = sum(expf, 1);
  
  arma::mat p(n, n_sources);
  
  
  
  for(int i=0; i<n; i++){
    for(int j =0; j<n_sources; j++){
      p(i,j) = expf(i,j)/sumexpf(i);
    }
  }
  
  return p;
  
}

//Log of likelihood added to prior
//
//
//[[Rcpp::export]]
double hcpp_pxr(int n_sources, int n_isotopes, int n_covariates,
            arma::vec d_prior,
            arma::mat x_scaled,
            arma::mat concentrationmeans, arma::mat sourcemeans,
            arma::mat correctionmeans,
            arma::mat corrsds, arma::mat sourcesds,
            arma::vec theta, arma::mat y, arma::vec c_prior,
            arma::mat sd_prior,
            arma::mat mu_prior,
            int uni_a_prior,
            int uni_b_prior){
  
  int n = y.n_rows;
  
  arma::mat beta(n_covariates, n_sources);
  
  for(int i = 0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
  }
  
  
  
  arma::mat p = hfn_pxr(theta, n_sources, n, n_covariates, x_scaled);
  
  // Initialize variables
  arma::mat mutotal(n, n_isotopes);
  arma::mat sigtotal(n, n_isotopes);
  
  // Calculate mutotal and sigtotal
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n_isotopes; j++) {
      double mutop = 0.0;
      double mubtm = 0.0;
      double sigsqtop = 0.0;
      double sigsqbtm = 0.0;
      
      for (int k = 0; k < n_sources; k++) {
        mutop += p(i, k) * concentrationmeans(k, j) * (sourcemeans(k, j) + correctionmeans(k, j));
        mubtm += p(i, k) * concentrationmeans(k, j);
        
        sigsqtop += pow(p(i, k) * concentrationmeans(k, j), 2) * (pow(sourcesds(k, j), 2) + pow(corrsds(k, j), 2));
        sigsqbtm += pow(p(i, k) * concentrationmeans(k, j), 2);
      }
      
      mutotal(i, j) = mutop / mubtm;
      // This is now sigtop/sigbottom x omicron (epsilon from mixsiar) + sigma_j
      sigtotal(i, j) = sqrt(sigsqtop / sigsqbtm * pow(exp(theta(j + n_sources * n_covariates +n_isotopes)), 0.5) + exp(theta(j + n_sources * n_covariates)));
    }
  }
  
  // Calculate hold
  double hold = 0.0;
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n_isotopes; j++) {
      hold += -log(sigtotal(i, j)) - 0.5 * log(2 * M_PI) -
        0.5 * pow((y(i, j) - mutotal(i, j)), 2) / pow(sigtotal(i, j), 2);
    }
  }
  
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < n_isotopes; j++) {
  //     hold += - n * n_isotopes * log(sigtotal(i, j)) - 0.5 * log(2 * M_PI) -
  //       0.5 * pow((y(i, j) - mutotal(i, j)), 2) / pow(sigtotal(i, j), 2);
  //   }
  // }
  
  // Calculate priornorm
  // double priornorm = 0.0;
  // 
  // // for(int i = 0; i<(n_covariates*n_sources + n_isotopes); i++){
  // //   priornorm += -(n_covariates * n_sources + n_isotopes) * log(sd_prior(i)) - 0.5 * log(2 * M_PI) -
  // //     0.5 * pow((theta(i)), 2) /pow(sd_prior(i),2);
  // // }
  // 
  // for(int i = 0; i<(n_covariates*n_sources + n_isotopes); i++){
  //   priornorm += -log(sd_prior(i)) - 0.5 * log(2 * M_PI) -
  //     0.5 * pow((theta(i)), 2) /pow(sd_prior(i),2);
  // }
  
  double betanorm = 0.0;
  
  //Zero here is prior for mu
  // for (int i = 0; i < n_covariates; i++) {
  //   for (int j = 0; j < n_sources; j++) {
  //     betanorm += -n_covariates * n_sources * log(1) - 0.5 * log(2 * M_PI) -
  //       0.5 * pow((beta(i, j) - 0), 2) / pow(1, 2);
  //   }
  // }
  
  for (int i = 0; i < n_covariates; i++) {
    for (int j = 0; j < n_sources; j++) {
      betanorm += - log(sd_prior(i,j)) - 0.5 * log(2 * M_PI) -
        0.5 * pow(beta(i, j) - mu_prior(i,j), 2) / pow(sd_prior(i,j), 2);
    }
  }
  // for (int i = 0; i < n_covariates; i++) {
  //   for (int j = 0; j < n_sources; j++) {
  //     betanorm += - log(1) - 0.5 * log(2 * M_PI) -
  //       0.5 * pow((beta(i, j) - 0), 2) / pow(1, 2);
  //   }
  // }
  
  // Calculate gammaprior
  double gammaprior = 0.0;
  
  for (int i = 0; i < n_isotopes; i++) {
    gammaprior += c_prior(i) * log(d_prior(i)) - lgamma(c_prior(i)) +
      (c_prior(i) - 1) * log(sqrt(exp(theta(i + n_sources * n_covariates)))) -
      d_prior(i) * sqrt(exp(theta(i + n_sources * n_covariates)));
  }
  
  // for (int i = 0; i < n_isotopes; i++) {
  //   gammaprior += -sqrt(exp(theta(i + n_sources * n_covariates)));
  // }
  
  double uniformprior = 1/(uni_b_prior-uni_a_prior);
  
  
  // Calculate totx
  //double totx = hold + priornorm;//gammaprior + betanorm + hold;
  double totx = gammaprior + betanorm + hold;
  return totx;
}

// // // This is basically the same as sim_theta but its using updated lambdas
// // // instead of set prior values

//[[Rcpp::export]]
double  log_q_cpp_pxr(arma::vec theta, arma::vec lambda,
                 int n_sources, int n_tracers, int S, int n_cov){
  
  int ncnsnt = n_sources * n_cov + (n_tracers * 2); // *2 for omicron and sigma
  
  arma::vec mean = lambda.subvec(0, ncnsnt - 1);
  
  int mat_size = ((ncnsnt) * (ncnsnt+1))/2;
  
  
  
  arma::vec sig = lambda.subvec(ncnsnt, ncnsnt + mat_size -1);
  
  
  
  arma::mat chol_var(ncnsnt, ncnsnt, arma::fill::zeros);
  
  int count = 0;
  
  for(int i = 0; i<ncnsnt; i++){
    for(int j = 0; j<ncnsnt; j++){
      if (j <= i){
        count +=1;
        chol_var((i),(j)) = sig(count-1);
        
      }
      
    }
  }
  
  // arma::mat prec = chol_prec * chol_prec.t();
  
  
  // arma::mat theta_mat(theta.memptr(), 1, n_sources * n_cov, false, true);
  // arma::mat lambda_mat(lambda.memptr(), 1, n_sources * n_cov, false, true);
  // arma::mat betaminusmean = theta_mat - lambda_mat;
  //
  //
  arma::mat thetaminusmean(ncnsnt, 1);
  
  for(int i = 0; i<ncnsnt; i++){
    
    thetaminusmean(i,0) = theta(i) - lambda(i);
  }
  
  // arma::mat tcholprec = chol_prec.t();
  // arma::mat Z = betaminusmean * tcholprec;
  // arma::mat tZ = Z.t();
  //
  //
  // NumericMatrix ZtZmat(1,1);
  //
  // ZtZmat = matmult(Z, tZ);
  //
  // double ZtZ = ZtZmat(0,0);
  
  
  
  // Z <- (y - mu)%*%cholPrec
  // return(- (m/2) * log(2* pi) + log(det(cholPrec)) - 0.5 * Z%*%t(Z))
  //
  
  arma::mat Z = arma::inv(chol_var) * thetaminusmean; //* betaminusmean.t();
  arma::mat ZtZ = Z.t() * Z;
  
  double ZtZ_scalar = ZtZ(0, 0);
  //arma::mat Sigma = chol_var * chol_var.t();
  //arma::double logdetSigma = log(det(Sigma));
  
  //NumericVector sig(n_cov*n_sources);
  
  
  // prod_sig = log(proddiag(tcholprec));
  
  // This is just the log determinant - so the log of the product of the diagonal
  double prod_sig = log(arma::prod(chol_var.diag()));
  
  double thetanorm = 0;
  
  //thetanorm = -  (ncnsnt/2) * log(2 * M_PI) - prod_sig - 0.5 * ZtZ_scalar;
  thetanorm = -  (ncnsnt/2) * log(2 * M_PI) - prod_sig - 0.5 * ZtZ_scalar;
  
  //thetanorm = -ncnsnt/2 * log(2 * M_PI) - 0.5 * logdetSigma - 0.5 *(thetaminusmean).t() *(Sigma/thetaminusmean);
  
  // double gamman = 0;
  // for (int i=0; i <(n_tracers); i++){
  //   gamman += lambda(mat_size + ncns +i) *
  //     log(lambda(mat_size + ncns +i + n_tracers))  -
  //     log(tgamma(lambda(mat_size + ncns +i)))  +
  //     (lambda(mat_size + ncns +i) - 1) * log(theta(i+ncns)) -
  //     lambda(mat_size + ncns +i + n_tracers) * theta((i+ncns));
  // }
  
  //Don't really understand this
  // arma::uvec indices = arma::regspace<arma::uvec>(mat_size + n_cov * n_sources, mat_size + n_cov * n_sources + n_tracers - 1);
  //
  // arma::vec gamman_vec = lambda.subvec(indices) % log(lambda.subvec(indices + n_tracers)) -
  //   lgamma(lambda.subvec(indices)) +
  //   (lambda.subvec(indices) - 1) % log(theta.subvec(arma::uvec(indices))) -
  //   lambda.subvec(indices + n_tracers) % theta.subvec(arma::uvec(arma::regspace(n_sources * n_cov, n_sources * n_cov + n_tracers - 1)));
  //
  // double gamman = arma::accu(gamman_vec);
  //
  
  //double x = sum_p+ gamman;
  //double x = thetanorm +gamman;
  
  return (thetanorm);
  
}


// getting the difference of hcpp_pxr and log_q_cpp_pxr
// [[Rcpp::export]]
double h_lambdacpp_pxr(int n_sources, int n_isotopes,
                   arma::vec beta_prior,
                   int n_covariates,
                   int S,
                   arma::mat concentrationmeans, arma::mat sourcemeans,
                   arma::mat correctionmeans,
                   arma::mat corrsds, arma::mat sourcesds,
                   arma::vec theta, arma::mat y,
                   arma::vec lambda,
                   arma::mat x_scaled,
                   arma::vec c_0,
                   arma::mat sd_prior,
                   arma::mat mu_prior,
                   int uni_a_prior,
                   int uni_b_prior) {
  
  return hcpp_pxr(n_sources, n_isotopes, n_covariates, beta_prior, x_scaled, concentrationmeans, sourcemeans, correctionmeans,
              corrsds, sourcesds, theta, y, c_0, sd_prior, mu_prior, uni_a_prior, uni_b_prior) - log_q_cpp_pxr(theta, lambda, n_sources, n_isotopes, S, n_covariates);
}




//Derivatives with respect to theta
// [[Rcpp::export]]
arma::vec delta_h_lambda_cpp_pxr(int n_sources, int n_tracers,
                             arma::vec beta_prior,
                             int n_covariates,
                             int S,
                             arma::mat concentrationmeans, arma::mat sourcemeans,
                             arma::mat correctionmeans,
                             arma::mat corrsds, arma::mat sourcesds,
                             arma::vec theta, arma::mat y,
                             arma::vec lambda,
                             arma::mat x_scaled,
                             arma::vec c_0,
                             arma::mat sd_prior,
                             arma::mat mu_prior,
                             double eps,
                             int uni_a_prior,
                             int uni_b_prior) {
  
  int ncnsnt = n_sources * n_covariates + n_tracers;
  int mat_size = ((ncnsnt) * (ncnsnt+1))/2;
  // eps = 0.001;
  
  double k = theta.n_elem;
  arma::vec ans(k);
  arma::vec d(k);
  arma::vec thetaplusd(k);
  arma::vec thetaminusd(k);
  
  arma::vec sig = lambda.subvec(ncnsnt, ncnsnt + mat_size -1);
  
  
  
  arma::mat chol_var(ncnsnt, ncnsnt, arma::fill::zeros);
  
  int count = 0;
  
  for(int i = 0; i<ncnsnt; i++){
    for(int j = 0; j<ncnsnt; j++){
      if (j <= i){
        count +=1;
        chol_var((i),(j)) = sig(count-1);
        
      }
      
    }
  }
  
  
  for(int i = 0; i<k; i++){
    
    for (int j = 0; j<k; j++){
      d(j) = 0;
    }
    d(i) = eps;
    
    
    for (int j = 0; j<k; j++){
      thetaplusd(j) = theta(j) + d(j);
      thetaminusd(j) = theta(j) - d(j);
    }
    ans(i) = (hcpp_pxr(n_sources, n_tracers,
              n_covariates,
              beta_prior,
              x_scaled,
              concentrationmeans, sourcemeans,
              correctionmeans,
              corrsds, sourcesds,
              thetaplusd, y,
              c_0,
              sd_prior,
              mu_prior,
              uni_a_prior,
              uni_b_prior) -
                hcpp_pxr(n_sources, n_tracers,
                     n_covariates,
                     beta_prior,
                     x_scaled,
                     concentrationmeans, sourcemeans,
                     correctionmeans,
                     corrsds, sourcesds,
                     thetaminusd, y,
                     c_0,
                     sd_prior,
                     mu_prior,
                     uni_a_prior,
                     uni_b_prior))/(2 * eps);
  }
  
  arma::vec log_q_deriv = - arma::inv(chol_var).t() * arma::inv(chol_var)* (theta - lambda.subvec(0, ncnsnt-1));
  
  // arma::vec log_q_deriv = - (chol_var) * (chol_var.t())* (theta - lambda.subvec(0, ncnsnt-1));
  arma::vec final_ans =  ans - log_q_deriv;
  
  return  final_ans;
}



//Nabla LB is the mean of delta_h_lambda times kappa_transpose
// [[Rcpp::export]]
arma::vec nabla_LB_cpp_pxr(arma::vec lambda, arma::mat theta,
                       int n_sources, int n_tracers, arma::vec beta_prior,
                       int S, int n_covariates,
                       arma::mat x_scaled,
                       arma::mat concentrationmeans,
                       arma::mat sourcemeans,
                       arma::mat correctionmeans,
                       arma::mat corrsds, arma::mat sourcesds,
                       arma::mat y,
                       arma::vec c_0,
                       arma::mat sd_prior,
                       arma::mat mu_prior,
                       arma::mat kappa,
                       int uni_a_prior,
                       int uni_b_prior){
  
  //int thetanrow = theta.n_rows; //THIS IS JUST S OOPS
  int lambdalength = lambda.n_elem;
  int ncnsnt = n_sources * n_covariates + n_tracers;
  int mat_size = ncnsnt * (ncnsnt +1)/2;
  arma::mat kappa_transpose = kappa.t();
  arma::vec ans(lambdalength);
  
  
  arma::mat big_delta_h_lambda(S, ncnsnt, arma::fill::zeros);
  
  
  for(int i = 0; i <S; i++){
    big_delta_h_lambda.row(i) = delta_h_lambda_cpp_pxr(n_sources, n_tracers,
                           beta_prior, n_covariates, S, concentrationmeans, sourcemeans,
                           correctionmeans,
                           corrsds,
                           sourcesds,
                           theta.row(i).t(), y, lambda, x_scaled, c_0, sd_prior, mu_prior,
                           0.0001, uni_a_prior, uni_b_prior).t();
  }
  
  
  // This creates the big_delta_h_lambda which is S * ncnsnt (transpose to get ncnsnt *S)
  // We take each col in turn (so we have ncnsnt * 1)
  // Multiply by each col of kappaT(1 * ncnsnt)
  // Get the vech of this (so then its mat_size)
  //Repeat S times
  
  arma::mat big_delta_h_lambda_transpose = big_delta_h_lambda.t();
  //This is ll *S
  
  //First thing we want to do is just get the average of this over cols
  // This is the nabla_mu_LB_lambda
  for (int i = 0; i < ncnsnt; i++) {
    ans(i) = arma::mean(big_delta_h_lambda_transpose.row(i));
  }
  
  
  // Now we want to take each col of big_delta_h_lambda_transpose (lt,1)
  // Multiply by each row of kappatranspose (1 * lt)
  // Get the vech of this (same code as reverse chol_prec thingy)
  //Loop over S
  //Save in an (lt, lt, S) array
  // Then vech each (lt,lt) - so we have matsize * S
  // Then average of all the vechs
  //Then output
  //
  arma::cube myArray(ncnsnt, ncnsnt, S, arma::fill::zeros);
  
  for(int i=0; i<S; i++){
    myArray.slice(i) = big_delta_h_lambda_transpose.col(i) * kappa_transpose.row(i);
  }
  
  // // Now we want to vech these so we get matsize * S
  //
  arma::mat vecharray(mat_size, S, arma::fill::zeros);
  
  
  for(int s = 0; s<S; s++){
    int count2 = 0;
    for(int i = 0; i<ncnsnt; i++){
      for(int j = 0; j<ncnsnt; j++){
        if (j <= i){
          count2 +=1;
          vecharray(count2-1, s) = myArray(i,j,s);
          
        }
        
      }
    }
  }
  //
  //
  for(int i = 0; i<mat_size; i++){
    
    ans(i+ncnsnt) = mean(vecharray.row(i));
    
  }
  
  
  return ans;
}








// estimate of LB
// [[Rcpp::export]]
double LB_lambda_cpp_pxr(arma::mat theta, arma::vec lambda,
                     int n_sources, int n_isotopes,
                     arma::vec beta_prior,
                     int n_covariates,
                     arma::mat x_scaled,
                     arma::mat concentrationmeans,
                     arma::mat sourcemeans,
                     arma::mat correctionmeans,
                     arma::mat corrsds,
                     arma::mat sourcesds,
                     arma::mat y,
                     arma::vec c_0,
                     arma::mat sd_prior,
                     arma::mat mu_prior,
                     int uni_a_prior,
                     int uni_b_prior){
  
  int S = theta.n_rows;
  
  arma::vec hlambdaapply(S);
  
  for(int i = 0; i <S; i++){
    hlambdaapply(i) = h_lambdacpp_pxr(n_sources, n_isotopes, beta_prior,
                 n_covariates, S,
                 concentrationmeans, sourcemeans,
                 correctionmeans, corrsds, sourcesds,
                 theta.row(i).t(), y, lambda,
                 x_scaled, c_0, sd_prior, mu_prior, uni_a_prior, uni_b_prior);
  }
  
  double ans = mean(hlambdaapply);
  
  return ans;
  
  
}


// Actually putting it all together and running it
// [[Rcpp::export]]
List run_VB_cpp_pxr(arma::vec lambdastart,
                int n_sources,
                int n_tracers,
                int n_covariates,
                int n,
                arma::vec beta_prior,
                arma::mat concentrationmeans,
                arma::mat sourcemeans,
                arma::mat correctionmeans,
                arma::mat corrsds,
                arma::mat sourcesds,
                arma::mat y,
                arma::mat x_scaled,
                int S,
                int P,
                double beta_1,
                double beta_2,
                int tau,
                double eps_0,
                int t_W,
                arma::vec c_prior,
                bool solo,
                arma::mat sd_prior,
                arma::mat mu_prior,
                int uni_a_prior,
                int uni_b_prior
){
  
  
  
  int lsl = lambdastart.n_elem;
  int ncnsnt = n_sources * n_covariates + n_tracers;
  
  //// NEED TO GENERATE EPS HERE
  // CAlling it kappa instead
  arma::mat kappa(ncnsnt, S, arma::fill::zeros);
  
  arma::vec mean_kappa(S, arma::fill::zeros);
  arma::mat var_kappa = arma::mat(S,S,arma::fill::eye);
  kappa = rMVNormCpp_pxr(ncnsnt, mean_kappa, var_kappa);
  
  //
  arma::mat theta = sim_thetacpp_pxr(S, lambdastart, n_sources, n_tracers, n_covariates, solo, kappa);
  //
  //
  //
  arma::vec g_0 = nabla_LB_cpp_pxr(lambdastart, theta,
                               n_sources, n_tracers,
                               beta_prior,
                               S, n_covariates, x_scaled,
                               concentrationmeans, sourcemeans,
                               correctionmeans, corrsds,
                               sourcesds, y, c_prior, sd_prior, mu_prior, kappa, uni_a_prior, uni_b_prior);
  
  
  
  arma::vec nu_0 = arma::square(g_0);
  
  arma::vec g_bar = g_0;
  
  arma::vec nu_bar = nu_0;
  
  arma::vec g_t(lsl);
  arma::vec nu_t(lsl);
  
  double patience = 0;
  bool stop = FALSE;
  double max_LB_bar = -arma::datum::inf;
  double alpha_t = 0;
  double t = 0;
  
  arma::vec LB(t_W + 1);
  LB.fill(NA_REAL); //(-std::numeric_limits<double>::infinity());
  
  arma::vec lambda = lambdastart;
  
  arma::vec VB_save(10000);
  arma::mat lambda_save(10000, lambda.n_elem);
  arma::vec h_lambda_cpp_save(10000);
  arma::vec lambda_max_LB(lambda.n_elem);
  //
  //
  while(!stop){
    
    //// NEED TO GENERATE EPS HERE
    // CAlling it kappa instead
    // arma::mat kappa(ncnsnt, S, arma::fill::zeros);
    //
    // arma::vec mean_kappa(S, arma::fill::zeros);
    // arma::mat var_kappa = arma::mat(S,S,arma::fill::eye);
    arma::mat kappa_loop(ncnsnt, S, arma::fill::zeros);
    kappa_loop = rMVNormCpp_pxr(ncnsnt, mean_kappa, var_kappa);
    
    theta = sim_thetacpp_pxr(S, lambda, n_sources, n_tracers, n_covariates, solo, kappa_loop);
    
    
    g_t = nabla_LB_cpp_pxr(lambda, theta, n_sources, n_tracers, beta_prior,
                       S, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds, sourcesds,
                       y, c_prior, sd_prior, mu_prior, kappa_loop, uni_a_prior, uni_b_prior);
    //
    //     // This is where we would have gradient clipping
    //     // So it would be some threshold/norm(g_t) * g_t
    //
    double norm_g_t = arma::norm(g_t, 2);
    if(norm_g_t >100){
      g_t = 100/norm_g_t * g_t;
    }
    //
    //
    //
    nu_t = arma::square(g_t);
    //
    //
    g_bar = beta_1 * g_bar + (1 - beta_1) * g_t;
    nu_bar = beta_2 * nu_bar + (1 - beta_2) * nu_t;
    //
    //
    arma::vec alpha_min(2);
    alpha_min(0) = eps_0;
    alpha_min(1) = eps_0 * (tau / (t + 1));
    //
    alpha_t = arma::min(alpha_min);
    //
    //     // h_lambda_cpp_save(t) = h_lambdacpp_pxr(n_sources, n_tracers, beta_prior,
    //     //                   n_covariates, S,
    //     //                   concentrationmeans, sourcemeans,
    //     //                   correctionmeans,
    //     //                   corrsds,sourcesds, theta.t(), y,
    //     //                   lambda, x_scaled, c_0, sd_prior);
    //
    //
    // Update lambda
    for(int i = 0; i<lsl; i++){
      
      lambda(i) = lambda(i) +  alpha_t * (g_bar(i) / std::sqrt(nu_bar(i)));
      
    }
    
    for(int i=0; i<lsl; i++){
      lambda_save(t,i) = lambda(i);
    }
    
    
    
    
    Rcout << "Iteration : " << t << "\n";
    
    
    //////////// This was written by Ahmed
    
    int r = t;
    int inn = 0;
    while(1){
      inn++;
      r = r/10;
      if(r == 0) break;
    }
    
    for(int j = 0 ; j < (13+inn) ;j++){
      Rcout<<"\b";
    }
    
    ///////////////////
    //
    //Compute the moving average LB if out of warm-up
    if(t<=t_W){
      
      
      LB(t) = LB_lambda_cpp_pxr(theta, lambda, n_sources, n_tracers,
         beta_prior,
         n_covariates, x_scaled,
         concentrationmeans, sourcemeans,
         correctionmeans,
         corrsds,sourcesds, y, c_prior, sd_prior, mu_prior, uni_a_prior, uni_b_prior);
    }
    else{
      for (int i = 0; i<(t_W-1); i++){
        LB(i) = LB(i+1);
      }
      
      
      LB(t_W) = LB_lambda_cpp_pxr(theta, lambda, n_sources, n_tracers,
         beta_prior,
         n_covariates, x_scaled,
         concentrationmeans, sourcemeans,
         correctionmeans,
         corrsds,sourcesds, y, c_prior, sd_prior, mu_prior, uni_a_prior, uni_b_prior);
      
      double LB_bar = arma::mean(LB);
      
      arma::vec maxbar(2);
      maxbar(0) = max_LB_bar;
      maxbar(1) = LB_bar;
      
      max_LB_bar = arma::max(maxbar);
      
      VB_save(t) = LB_bar;
      
      
      
      
      
      if(LB_bar>= max_LB_bar){
        patience = 0;
        lambda_max_LB = lambda;
      } else{
        patience = patience +1;
      }
      
    }
    if(patience>P){
      stop = TRUE;
    }
    
    if(t>10000){
      stop = TRUE;
    }
    
    t = t + 1;
  }
  
  
  return Rcpp::List::create(Rcpp::Named("kappa") = kappa,
                            Rcpp::Named("theta") = theta,
                            Rcpp::Named("g_o") = g_0 ,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("mean_LB_bar") = VB_save,
                            Rcpp::Named("lambda_save") = lambda_save,
                            Rcpp::Named("iteration") = t,
                            Rcpp::Named("h_lambda") = h_lambda_cpp_save,
                            Rcpp::Named("t") = t,
                            Rcpp::Named("LB") = LB,
                            Rcpp::Named("g_t") = g_t,
                            Rcpp::Named("lambda_max_LB") = lambda_max_LB
  );
  
}


