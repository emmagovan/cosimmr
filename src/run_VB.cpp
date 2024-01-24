#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


// [[Rcpp::export]]
arma::vec digamma_wrapper(const arma::vec& x) {
  Rcpp::Function digamma_R = Rcpp::Function("digamma");
  return Rcpp::as<arma::vec>(digamma_R(x));
}



// [[Rcpp::export]]
arma::mat rMVNormCpp(int n, arma::vec Mean, arma::mat Var) {
  int ncols = Var.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(Mean, 1, n).t() + Y * arma::chol(Var);
}


//simulates S samples of theta - so several rnorms + rgamma
//[[Rcpp::export]]
arma::mat sim_thetacpp(int S, arma::vec lambda, int n_sources,
                       int n_tracers, int n_cov, bool solo){
  int ncns = n_sources * n_cov;
  
  arma::mat theta(S, (ncns + n_tracers));
  
  arma::vec mean_beta = lambda.subvec(0, ncns -1);
  
  
  int mat_size = ((ncns) * (ncns+1))/2;
  
  arma::vec sig_beta = lambda.subvec(ncns, ncns + mat_size -1);
  
  arma::mat chol_prec(ncns, ncns, arma::fill::zeros);
  
  int count = 0;
  for(int i = 0; i<ncns; i++){
    for(int j = 0; j<ncns; j++){
      if (j <= i){
        count +=1;
        chol_prec((i),(j)) = sig_beta(count-1);
      }
    }
  }
  
  arma::mat normmat(S, ncns);
  arma::mat prec = chol_prec * chol_prec.t();
  // arma::mat solve_prec = arma::inv(prec);
  // arma::mat var = solve_prec.t();
  arma::mat var = arma::inv(prec);
  
  normmat = rMVNormCpp(S, mean_beta, var);
  
  for(int i=0; i<ncns; i++){
    theta.col(i) = normmat.col(i);
  }
  
  
  arma::vec solovec(S, arma::fill::ones);
  solovec *= 1000;
  
  if(!solo){
    for(int i = 0; i<n_tracers; i++){
      // shape and scale are params here
      theta.col(i+ncns) = as<arma::vec>(Rcpp::rgamma(S,  lambda(mat_size + ncns +i),
                                        //lambda(mat_size + ncns +i + n_tracers)));
                                        1/lambda(mat_size + ncns +i + n_tracers)));
    }
  } else {
    for(int i = 0; i<n_tracers; i++){
      theta.col(i+ncns) = solovec + 0.00001 * as<arma::vec>(Rcpp::rgamma(S,  lambda(mat_size + ncns +i),
                                                            1/lambda(mat_size + ncns +i + n_tracers)));
    }
  }
  
  
  return theta;
}



// calculates p - which is exp(f) / sum exp(f)
// and f is X * B
// X is a matrix and B is a vector

//[[Rcpp::export]]
arma::mat hfn(arma::vec theta,
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
// //
// //
//Log of likelihood added to prior
//
//
//[[Rcpp::export]]
double hcpp(int n_sources, int n_isotopes, int n_covariates,
            arma::vec d_prior,
            arma::mat x_scaled,
            arma::mat concentrationmeans, arma::mat sourcemeans,
            arma::mat correctionmeans,
            arma::mat corrsds, arma::mat sourcesds,
            arma::vec theta, arma::mat y, arma::vec c_prior,
            arma::mat sd_prior){
  
  int n = y.n_rows;
  
  arma::mat beta(n_covariates, n_sources);
  
  for(int i = 0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
  }
  
  
  
  arma::mat p = hfn(theta, n_sources, n, n_covariates, x_scaled);
  
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
      sigtotal(i, j) = sqrt((sigsqtop / sigsqbtm + 1/theta(j + n_sources * n_covariates)));
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
  
  // Calculate betanorm
  double betanorm = 0.0;
  
  //Zero here is prior for mu
  for (int i = 0; i < n_covariates; i++) {
    for (int j = 0; j < n_sources; j++) {
      betanorm += -n_covariates * n_sources * log(sd_prior(i,j)) - 0.5 * log(2 * M_PI) -
        0.5 * pow((beta(i, j) - 0), 2) / pow(sd_prior(i,j), 2);
    }
  }
  
  // for (int i = 0; i < n_covariates; i++) {
  //   for (int j = 0; j < n_sources; j++) {
  //     betanorm += - log(sd_prior(i,j)) - 0.5 * log(2 * M_PI) -
  //       0.5 * pow((beta(i, j) - 0), 2) / pow(sd_prior(i,j), 2);
  //   }
  // }
  
  // Calculate gammaprior
  double gammaprior = 0.0;
  
  for (int i = 0; i < n_isotopes; i++) {
    gammaprior += c_prior(i) * log(d_prior(i)) - lgamma(c_prior(i)) +
      (c_prior(i) - 1) * log(theta(i + n_sources * n_covariates)) -
      d_prior(i) * theta(i + n_sources * n_covariates);
  }
  
  // Calculate totx
  double totx = gammaprior + betanorm + hold;
  
  return totx;
}

// // // This is basically the same as sim_theta but its using updated lambdas
// // // instead of set prior values

//[[Rcpp::export]]
double log_q_cpp(arma::vec theta, arma::vec lambda,
                 int n_sources, int n_tracers, int S, int n_cov){
  
  int ncns = n_sources * n_cov;
  
  arma::vec mean_beta = lambda.subvec(0, ncns - 1);
  
  int mat_size = ((ncns) * (ncns+1))/2;
  
  
  
  arma::vec sig_beta = lambda.subvec(ncns, ncns + mat_size -1);
  
  
  
  arma::mat chol_prec(ncns, ncns, arma::fill::zeros);
  
  int count = 0;
  
  for(int i = 0; i<ncns; i++){
    for(int j = 0; j<ncns; j++){
      if (j <= i){
        count +=1;
        chol_prec((i),(j)) = sig_beta(count-1);
        
      }
      
    }
  }
  
  arma::mat prec = chol_prec * chol_prec.t();
  
  
  // arma::mat theta_mat(theta.memptr(), 1, n_sources * n_cov, false, true);
  // arma::mat lambda_mat(lambda.memptr(), 1, n_sources * n_cov, false, true);
  // arma::mat betaminusmean = theta_mat - lambda_mat;
  //
  //
  arma::mat betaminusmean(1, ncns);
  
  for(int i = 0; i<ncns; i++){
    
    betaminusmean(0,i) = theta(i) - lambda(i);
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
  
  arma::mat Z = betaminusmean * chol_prec; //* betaminusmean.t();
  arma::mat ZtZ = Z * Z.t();
  
  double ZtZ_scalar = ZtZ(0, 0);
  
  
  //NumericVector sig(n_cov*n_sources);
  
  
  // prod_sig = log(proddiag(tcholprec));
  double prod_sig = log(arma::prod(chol_prec.diag()));
  
  double thetanorm = 0;
  
  thetanorm = -  (ncns/2) * log(2 * M_PI) + prod_sig - 0.5 * ZtZ_scalar;
  
  
  double gamman = 0;
  for (int i=0; i <(n_tracers); i++){
    gamman += lambda(mat_size + ncns +i) *
      log(lambda(mat_size + ncns +i + n_tracers))  -
      log(tgamma(lambda(mat_size + ncns +i)))  +
      (lambda(mat_size + ncns +i) - 1) * log(theta(i+ncns)) -
      lambda(mat_size + ncns +i + n_tracers) * theta((i+ncns));
  }
  
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
  double x = thetanorm +gamman;
  
  return (x);
  
}

// [[Rcpp::export]]
arma::vec delta_lqltcppauto(arma::vec lambda, arma::vec theta,
                            int n_sources,
                            int n_tracers, int n_covariates,
                            int S){
  
  int k = lambda.n_elem;
  
  arma::vec ans(k, arma::fill::zeros);
  
  int ncns = n_sources * n_covariates;
  
  arma::vec mean_beta = lambda.subvec(0, ncns - 1);
  
  int mat_size = ((ncns) * (ncns+1))/2;
  
  arma::vec sig_beta = lambda.subvec(ncns, ncns + mat_size -1);
  
  arma::mat chol_prec(ncns, ncns, arma::fill::zeros);
  
  int count = 0;
  
  for(int i = 0; i<ncns; i++){
    for(int j = 0; j<ncns; j++){
      if (j <= i){
        count +=1;
        chol_prec((i),(j)) = sig_beta(count-1);
        
      }
      
    }
  }
  
  arma::mat prec = chol_prec * chol_prec.t();
  arma::mat var = arma::inv(prec);
  
  // This is the mean part
  ans.subvec(0, ncns - 1) = prec * (theta.subvec(0, ncns -1) - mean_beta);
  
  //// This is a symmetric matrix - just want upper triangular elements and want to format as a vector
  // arma::mat ans2 = (-0.5 *(prec - prec * (theta.subvec(0, ncns -1) - mean_beta) * (theta.subvec(0, ncns -1) - mean_beta).t() * prec)); //* chol_prec;
  arma::mat ans2(ncns, ncns);
  // For the sigma part - need to fill different parts
  // Can do the matrix and then the off diagonals
  
  // Derivatives from here
  // https://www.sciencedirect.com/science/article/pii/S2452306222000168#sec0032
  //chol_prec = L^-1^T
  //prec = L^-1TL^-1
  
  arma::vec ytilde = (theta.subvec(0, ncns -1) - mean_beta);
  // double count3 = 0;
  //arma::vec test(4);
  for(int i = 0; i<ncns; i++){
    for(int j = 0; j<ncns; j++){
      if(i == j){
        double sumylambda = 0;
        for(int m=0; m<ncns; m++){
          
          sumylambda += theta(m) * chol_prec(m,i);
        }
        ans2(i,j) = 1/chol_prec(i,j) -ytilde(i) * sumylambda;
      }
    }}
  
  for(int i = 0; i<ncns; i++){
    for(int j = 0; j<ncns; j++){
      if(j < i){
        double sumylambda = 0;
        for(int m=0; m<ncns; m++){
          
          sumylambda += theta(m) * chol_prec(m,j);
        }
        ans2(i,j) = -ytilde(i) * sumylambda;
      }
    }}
  
  arma::vec derivative_sigma(mat_size);
  
  int count2 = 0;
  //This is just to extract the upper triangular elements in the right order
  
  for(int i = 0; i<ncns; i++){
    for(int j = 0; j<ncns; j++){
      if (j <= i){
        count2 +=1;
        derivative_sigma(count2-1) = ans2(i,j);
        
      }
      
    }
  }
  
  ans.subvec(ncns, ncns + mat_size -1) = derivative_sigma;
  
  arma::vec c = lambda.subvec(ncns + mat_size, ncns + mat_size + n_tracers - 1);
  arma::vec d = lambda.subvec(ncns + mat_size + n_tracers, ncns + mat_size + 2 * n_tracers-1);
  
  // When I did this before I didn't need the /lgamma - not sure why, normally included in the digamma for some reason
  //ans.subvec(ncns + mat_size, ncns + mat_size + n_tracers - 1) = log(d) - digamma_wrapper(c)/lgamma(c) - log(theta(ncns));
  ans.subvec(ncns + mat_size, ncns + mat_size + n_tracers - 1) = log(d) - digamma_wrapper(c) + log(theta.subvec(ncns, ncns+1));
  
  
  ans.subvec(ncns + mat_size+ n_tracers, ncns + mat_size + 2 * n_tracers - 1) = c/d - theta.subvec(ncns, ncns+1);
  
  
  return ans;
  
}


// getting the difference of hcpp and log_q_cpp
// [[Rcpp::export]]
double h_lambdacpp(int n_sources, int n_isotopes,
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
                   arma::mat sd_prior) {
  
  return hcpp(n_sources, n_isotopes, n_covariates, beta_prior, x_scaled, concentrationmeans, sourcemeans, correctionmeans,
              corrsds, sourcesds, theta, y, c_0, sd_prior) - log_q_cpp(theta, lambda, n_sources, n_isotopes, S, n_covariates);
}

//calculating covariance of 2 matrices
// [[Rcpp::export]]
arma::mat cov_mat_cpp(arma::mat x_arma, arma::mat y_arma) {
  
  int xrow = x_arma.n_rows;
  
  
  arma::rowvec meanx = arma::mean(x_arma, 0);
  arma::rowvec meany = arma::mean(y_arma, 0);
  
  arma::mat xminusmean = x_arma.each_row() - meanx;
  arma::mat yminusmean = y_arma.each_row() - meany;
  
  arma::mat covmat = (xminusmean.t() * yminusmean) / (xrow - 1);
  
  return covmat;
}



//Nabla LB is the mean of delta_lqlt element-wise multiplied by h_lambda
// [[Rcpp::export]]
arma::colvec nabla_LB_cpp(arma::vec lambda, arma::mat theta,
                          int n_sources, int n_tracers, arma::vec beta_prior,
                          int S, int n_covariates,
                          arma::mat x_scaled,
                          arma::mat concentrationmeans,
                          arma::mat sourcemeans,
                          arma::mat correctionmeans,
                          arma::mat corrsds, arma::mat sourcesds,
                          arma::mat y,
                          arma::vec c,
                          arma::vec c_0,
                          arma::mat sd_prior){
  
  int thetanrow = theta.n_rows;
  int lambdalength = lambda.n_elem;
  
  arma::mat big_c(thetanrow, c.n_elem);
  
  arma::mat big_delta_lqlt(thetanrow, lambdalength, arma::fill::zeros);
  arma::vec big_h_lambda(thetanrow);
  arma::mat big_h_lambda_rep(lambdalength, thetanrow);
  arma::mat big_h_lambda_rep_transpose(thetanrow, lambdalength);
  
  
  for(int i = 0; i <thetanrow; i++){
    // big_delta_lqlt.row(i) = delta_lqltcpp(lambda, theta.row(i).t(), 0.001, n_sources, n_tracers,
    //                   n_covariates, S).t();
    big_delta_lqlt.row(i) = delta_lqltcppauto(lambda, theta.row(i).t(), n_sources, n_tracers,
                       n_covariates, S).t();
  }
  
  for(int i =0; i<thetanrow; i++){
    big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers, beta_prior,
                 n_covariates, S,
                 concentrationmeans, sourcemeans,
                 correctionmeans,
                 corrsds,sourcesds, theta.row(i).t(), y,
                 lambda, x_scaled, c_0, sd_prior);
  }
  
  
  for (int i = 0; i < lambdalength; i++) {
    big_h_lambda_rep.row(i) = big_h_lambda.t();
  }
  
  big_h_lambda_rep_transpose = big_h_lambda_rep.t();
  
  for (int i = 0; i < thetanrow; i++) {
    big_c.row(i) = c.t();
  }
  
  arma::mat big_h_minus_c = big_h_lambda_rep_transpose - big_c;
  
  arma::mat ansmat(big_delta_lqlt.n_rows, big_h_minus_c.n_cols);
  
  for (int i = 0; i < big_delta_lqlt.n_rows; i++) {
    for (int j = 0; j < big_delta_lqlt.n_cols; j++) {
      ansmat(i, j) = big_delta_lqlt(i, j) * big_h_minus_c(i, j);
    }
  }
  
  arma::vec ans(lambdalength);
  
  for (int i = 0; i < ansmat.n_cols; i++) {
    ans(i) = arma::mean(ansmat.col(i));
  }
  
  
  
  
  return ans;
}


// calculate control variate (big formula)
// [[Rcpp::export]]
arma::vec control_var_cpp(arma::vec lambda,
                          arma::mat theta,
                          int n_sources, int n_tracers,
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
                          arma::mat sd_prior){
  
  int S = theta.n_rows;
  
  int lambdallength = lambda.n_elem;
  
  
  arma::mat big_delta_lqlt(S, lambdallength);
  arma::mat big_h_lambda_rep(lambdallength, S);
  arma::mat big_h_lambda_rep_transpose(S, lambdallength);
  arma::vec big_h_lambda(S);
  arma::vec big_h_lambda_transpose(S);
  
  for(int i = 0; i <S; i++){
    //big_delta_lqlt.row(i) = delta_lqltcpp(lambda, theta.row(i).t(), 0.001, n_sources, n_tracers,
    //                   n_covariates, S).t();
    big_delta_lqlt.row(i) = delta_lqltcppauto(lambda, theta.row(i).t(), n_sources, n_tracers,
                       n_covariates, S).t();
  }
  
  for(int i =0; i<S; i++){
    big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers, beta_prior,
                 n_covariates, S,
                 concentrationmeans, sourcemeans,
                 correctionmeans,
                 corrsds,sourcesds, theta.row(i).t(), y,
                 lambda, x_scaled, c_0, sd_prior);
  }
  
  
  for(int i =0; i<lambdallength; i++){
    big_h_lambda_rep.row(i) = big_h_lambda.t();
  }
  
  
  big_h_lambda_rep_transpose = big_h_lambda_rep.t();
  
  arma::mat big_nabla(S, lambdallength);
  
  for (int i = 0; i < S; i++)
  {
    for (int j = 0; j < lambdallength; j++) {
      
      
      big_nabla(i,j) = big_delta_lqlt(i,j) * big_h_lambda_rep_transpose(i,j);
      
      
    }
  }
  
  arma::vec var_big_delta_lqlt(lambdallength);
  
  for(int i = 0; i<lambdallength; i++){
    var_big_delta_lqlt(i) = arma::var(big_delta_lqlt.col(i));
  }
  
  arma::mat covmat = cov_mat_cpp(big_nabla, big_delta_lqlt);
  
  arma::vec diag(lambdallength);
  for(int i =0; i<lambdallength; i++){
    for(int j =0; j<lambdallength; j++){
      if(i == j){
        diag(i) = covmat(i,j);
      }
    }}
  
  arma::vec ans(lambdallength);
  for(int i =0; i<lambdallength; i++){
    ans(i) = diag(i)/var_big_delta_lqlt(i);
  }
  
  return ans;
}

// estimate of LB
// [[Rcpp::export]]
double LB_lambda_cpp(arma::mat theta, arma::vec lambda,
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
                     arma::mat sd_prior){
  
  int S = theta.n_rows;
  
  arma::vec hlambdaapply(S);
  
  for(int i = 0; i <S; i++){
    hlambdaapply(i) = h_lambdacpp(n_sources, n_isotopes, beta_prior,
                 n_covariates, S,
                 concentrationmeans, sourcemeans,
                 correctionmeans, corrsds, sourcesds,
                 theta.row(i).t(), y, lambda,
                 x_scaled, c_0, sd_prior);
  }
  
  double ans = mean(hlambdaapply);
  
  return ans;
  
  
}


// Actually putting it all together and running it
// [[Rcpp::export]]
List run_VB_cpp(arma::vec lambdastart,
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
                arma::mat sd_prior
){
  
  
  
  int lsl = lambdastart.n_elem;
  
  arma::mat theta = sim_thetacpp(S, lambdastart, n_sources, n_tracers, n_covariates, solo);
  
  arma::vec c = control_var_cpp(lambdastart, theta, n_sources, n_tracers, beta_prior,
                                n_covariates, x_scaled,
                                concentrationmeans,
                                sourcemeans, correctionmeans,
                                corrsds, sourcesds, y, c_prior, sd_prior);
  
  arma::vec c_0(lsl, arma::fill::zeros);
  
  arma::vec g_0 = nabla_LB_cpp(lambdastart, theta,
                               n_sources, n_tracers,
                               beta_prior,
                               S, n_covariates, x_scaled,
                               concentrationmeans, sourcemeans,
                               correctionmeans, corrsds,
                               sourcesds, y, c_0, c_prior, sd_prior);
  
  
  
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
  
  
  while(!stop){
    
    theta = sim_thetacpp(S, lambda, n_sources, n_tracers, n_covariates, solo);
    
    
    g_t = nabla_LB_cpp(lambda, theta, n_sources, n_tracers, beta_prior,
                       S, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds, sourcesds,
                       y, c, c_prior, sd_prior);
    
    // This is where we would have gradient clipping
    // So it would be some threshold/norm(g_t) * g_t
    
    // double norm_g_t = arma::norm(g_t, 2);
    // if(norm_g_t >100){
    //   g_t = 100/norm_g_t * g_t;
    // }
    
    
    c = control_var_cpp(lambda, theta,n_sources,n_tracers, beta_prior,
                        n_covariates, x_scaled,
                        concentrationmeans, sourcemeans,
                        correctionmeans,
                        corrsds,sourcesds, y, c_prior, sd_prior);
    
    
    nu_t = arma::square(g_t);
    
    
    g_bar = beta_1 * g_bar + (1 - beta_1) * g_t;
    nu_bar = beta_2 * nu_bar + (1 - beta_2) * nu_t;
    
    
    arma::vec alpha_min(2);
    alpha_min(0) = eps_0;
    alpha_min(1) = eps_0 * (tau / (t + 1));
    
    alpha_t = arma::min(alpha_min);
    
    
    // Update lambda
    for(int i = 0; i<lsl; i++){
      
      lambda(i) = lambda(i) +  alpha_t * (g_bar(i) / std::sqrt(nu_bar(i)));
      
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
    
    //Compute the moving average LB if out of warm-up
    if(t<=t_W){
      
      
      LB(t) = LB_lambda_cpp(theta, lambda, n_sources, n_tracers,
         beta_prior,
         n_covariates, x_scaled,
         concentrationmeans, sourcemeans,
         correctionmeans,
         corrsds,sourcesds, y, c_prior, sd_prior);
    }
    else{
      for (int i = 0; i<(t_W-1); i++){
        LB(i) = LB(i+1);
      }
      
      
      LB(t_W) = LB_lambda_cpp(theta, lambda, n_sources, n_tracers,
         beta_prior,
         n_covariates, x_scaled,
         concentrationmeans, sourcemeans,
         correctionmeans,
         corrsds,sourcesds, y, c_prior, sd_prior);
      
      double LB_bar = arma::mean(LB);
      
      arma::vec maxbar(2);
      maxbar(0) = max_LB_bar;
      maxbar(1) = LB_bar;
      
      max_LB_bar = arma::max(maxbar);
      
      VB_save(t) = LB_bar;
      
      for(int i=0; i<lambda.n_elem; i++){
        lambda_save(t,i) = lambda(i);
      }
      
      if(LB_bar>= max_LB_bar){
        patience = 0;
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
  
  
  return Rcpp::List::create(Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("mean_LB") = VB_save,
                            Rcpp::Named("lambda_save") = lambda_save,
                            Rcpp::Named("iteration") = t
  );
  
}
