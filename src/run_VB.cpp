#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

static double const log2pi = std::log(2.0 * M_PI);


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// this is calculating the determinant specifically of an upper triangular matrix
// via calculating the product of the diagonal
// [[Rcpp::export]]
double proddiag(NumericMatrix x) { 
  double a = x.ncol();
  double b = x.nrow();
  double ans = 0;
  for(int i=0; i<a; i++){
    for(int j = 0; j<b; j++){
      if(i==j){
        ans += log(x(i,j));
        
      }
    }
  }
  return(exp(ans));
}



// 
// //dnorm 
// // [[Rcpp::export]]
// arma::vec dmvnrm_arma_fast(arma::mat const &x,  
//                            arma::rowvec const &mean,  
//                            arma::mat const &sigma, 
//                            bool const logd = true) { 
//   using arma::uword;
//   uword const n = x.n_rows, 
//     xdim = x.n_cols;
//   arma::vec out(n);
//   arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
//   double const rootisum = arma::sum(log(rooti.diag())), 
//     constants = -(double)xdim/2.0 * log2pi, 
//     other_terms = rootisum + constants;
//   
//   arma::rowvec z;
//   for (uword i = 0; i < n; i++) {
//     z = (x.row(i) - mean);
//     inplace_tri_mat_mult(z, rooti);
//     out(i) = other_terms - 0.5 * arma::dot(z, z);     
//   }  
//   
//   if (logd)
//     return out;
//   return exp(out);
// }


//matrix multiplication
// [[Rcpp::export]]
NumericMatrix matmult(NumericMatrix x, NumericMatrix y) {
  NumericMatrix ans(x.nrow(), y.ncol());
  
  for(int i = 0; i<x.nrow(); i++){
    for(int j=0; j<y.ncol(); j++){
      for(int k =0; k<y.nrow(); k++){
        
        ans(i,j) += x(i,k) * y(k,j);
      }
    }}
  
  return ans;
}


//for our specific q we are technically using the transpose of the chol...
// [[Rcpp::export]]
NumericMatrix crossprod(NumericMatrix X){
  
  NumericMatrix tX(X.ncol(), X.nrow());
  
  for(int i=0; i<X.nrow(); i++){
    for (int j=0; j < X.ncol(); j++){
      tX(j,i) = X(i,j);
    }}
  
  
  NumericMatrix ans(X.nrow(), X.ncol());
  
  for(int i = 0; i<X.ncol(); i++){
    for(int j=0; j<X.ncol(); j++){
      for(int n =0; n<X.nrow(); n++){
        
        ans(i,j) += tX(i,n) * X(n,j);
      }
    }}
  return(ans);
}


//rnorm sample


// [[Rcpp::export]]
NumericMatrix rMVNormCpp(int n, arma::vec Mean, arma::mat Var) {
  int ncols = Var.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return Rcpp::wrap(arma::repmat(Mean, 1, n).t() + Y * arma::chol(Var));
}

// 
// // [[Rcpp::export]]
// NumericMatrix rMVNormCpp(const double n,
//                          const arma::vec mu,
//                          const NumericMatrix U) {
//   
//   
//   // Dimension of MVN
//   int p = mu.size();
//   
//   // Simulate iid standard normals
//   arma::mat Z(p, n);
//   Z.imbue(norm_rand);
//   
//   // Now backsolve and add back on the means
//   arma::mat X = solve(as<arma::mat>(U), Z);
//   for ( int i = 0; i < n; ++i ) {
//     X.col(i) += mu;
//   }
//   
//   return Rcpp::wrap(X.t());
// }


//solves matrix - returns transpose though
// [[Rcpp::export]]
NumericMatrix solvearma(const NumericMatrix X) {
  
  arma::mat b = arma::eye(X.nrow(), X.ncol());
  
  
  // Now backsolve and add back on the means
  arma::mat ans = solve(as<arma::mat>(X), b);
  
  
  return Rcpp::wrap(ans.t());
}

//simulates S samples of theta - so several rnorms + rgamma
//[[Rcpp::export]]
NumericMatrix sim_thetacpp(int S, NumericVector lambda, int n_sources,
                           int n_tracers, int n_cov, bool solo){
  
  NumericMatrix theta(S, (n_cov*n_sources + n_tracers));
  
  NumericVector mean_beta(n_cov*n_sources);
  
  int mat_size = ((n_cov*n_sources) * (n_cov*n_sources+1))/2;
  
  
  for(int i=0; i<n_cov*n_sources; i++){
    mean_beta(i) =  lambda(i);
  }
  
  
  
  NumericVector sig_beta(mat_size);
  
  
  
  for(int i = 0; i<mat_size; i++){
    sig_beta(i) = lambda(i + n_cov*n_sources);
    
  }
  
  
  NumericMatrix chol_prec(n_cov * n_sources, n_cov * n_sources);
  
  int count = 0;
  for(int i = 0; i<n_sources* n_cov; i++){
    for(int j = 0; j<n_sources * n_cov; j++){
      if (i <= j){
        count +=1;
        chol_prec((i),(j)) = sig_beta(count-1);
        
        
      }
      else{
        chol_prec(i,j) = 0;
      }
      
    }
  }
  
  NumericMatrix normmat(S, n_sources*n_cov);
  NumericMatrix prec(n_sources*n_cov, n_sources*n_cov);
  prec = crossprod(chol_prec);
  
  NumericMatrix solve_prec(n_sources*n_cov, n_sources*n_cov);
  solve_prec = solvearma(prec);
  
  // solvearma returns the transpose of the inverse of the precision matrix
  //so need to transpose it back
  NumericMatrix var(n_sources*n_cov, n_sources*n_cov);
  for(int i=0; i<n_sources*n_cov; i++){
    for(int j=0; j<n_sources*n_cov; j++){
      var(i,j) = solve_prec(j,i);
    }
  }
  
  normmat = rMVNormCpp(S, mean_beta, as<arma::mat>(var));
  
  for(int i=0; i<n_sources*n_cov; i++){
    theta(_,i) = normmat(_,i);
  }
  
  // arma::mat theta_arma(S, (n_cov*n_sources + n_tracers)); //Want to go from chol_prec array to
  // // A matrix of thetas generated using rMVNormCpp
  // 
  // for(int i=0; i<n_cov; i++){
  //   theta_arma.submat(0, (i)*n_sources, S-1, (i+1)*n_sources - 1) = as<arma::mat>(rMVNormCpp(S, mean_beta_mat(i,_), Rcpp::wrap(chol_prec.slice(i))));
  // }
  // 
  // theta = Rcpp::wrap(theta_arma);
  // 
  NumericVector solovec(S);
  for(int s = 0; s<S; s++){
    solovec(s) = 1000;
  }
  
  if(solo == FALSE) for(int i = 0; i<n_tracers; i++){
    theta(_,i+n_sources*n_cov) = (Rcpp::rgamma(S,  lambda(mat_size + n_cov *n_sources +i),
                                  1/lambda(mat_size + n_cov *n_sources +i + n_tracers)));
  }
  else if(solo == TRUE)for(int i = 0; i<n_tracers; i++){
    theta(_,i+n_sources*n_cov) = solovec + 0.00001 * (Rcpp::rgamma(S,  lambda(mat_size + n_cov *n_sources +i),
                                                      1/lambda(mat_size + n_cov *n_sources +i + n_tracers)));
  }
  
  
  return theta;
  
  // so this is K betas and then n_isotope taus 
  // Tau is gamma distributed with shape and rate
  // Need to check I'm using tau the whole way through and not sigma
}



//calculates p - which is exp(f) / sum exp(f)
//and f is X * B
//X is a matrix and B is a vector

//[[Rcpp::export]]
NumericMatrix hfn(NumericVector theta, 
                  int n_sources, int n, int n_cov, 
                  NumericMatrix x_scaled){
  NumericMatrix p(n, n_sources);
  NumericMatrix exptheta(n_sources, n_sources);
  NumericMatrix f(n, n_sources);
  NumericMatrix beta(n_cov, n_sources);
  
  
  for(int i = 0; i<n_cov; i++){
    for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
  }
  
  f = matmult(x_scaled, beta);
  
  
  NumericMatrix expf(n, n_sources); 
  
  for(int i =0; i<n; i++){
    for(int j = 0; j<n_sources; j++){
      expf(i,j) = exp(f(i,j));
    }
  }
  
  NumericVector sumexpf(n);
  
  for(int i = 0; i<n; i++){
    for(int j=0; j<n_sources; j++){
      sumexpf(i) +=expf(i,j);
    }
  }
  
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
double hcpp(int n_sources, int n_isotopes, int n_covariates,
            NumericVector beta_prior,
            NumericMatrix x_scaled,
            NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
            NumericMatrix correctionmeans,
            NumericMatrix corrsds, NumericMatrix sourcesds,
            NumericVector theta, NumericMatrix y, NumericVector c_0){
  
  
  
  Rcpp::NumericMatrix beta(n_covariates, n_sources);
  
  for(int i = 0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
  }
  
  int n = y.rows();
  
  NumericMatrix p(n, n_sources);
  
  p = hfn(theta, n_sources, n, n_covariates, x_scaled);
  
  
  
  // Setting prior values for hyper parameters
  NumericMatrix prior_means(n_covariates, n_sources);
  NumericMatrix prior_sd(n_covariates, n_sources);
  //At some stage change these to be consistent
  //NumericVector c_0(n_isotopes);
  NumericVector d_0(n_isotopes);
  
  // Setting up prior values
  for(int i=0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
      
      prior_means(i,j) = 0;
      prior_sd(i,j) = 1;
    }
  }
  
  
  for (int i = 0; i<n_isotopes; i++){
    //   //c_0(i) = 0.01;
    d_0(i) = beta_prior(i);
  }
  
  // Make sure these are set to zero
  NumericMatrix mutop(n,n_isotopes);
  NumericMatrix mubtm(n,n_isotopes);
  
  for(int j=0; j<n_isotopes; j++){
    for(int i=0; i<n; i++){
      for(int k=0; k<n_sources; k++){
        mutop(i,j) += p(i,k)*concentrationmeans(k,j) * (sourcemeans(k,j) + correctionmeans(k,j));
        mubtm(i,j) += p(i,k)*concentrationmeans(k,j);
      }
    }
  }
  
  
  NumericMatrix mutotal(n,n_isotopes);
  // 
  for(int i=0; i<(n); i++){
    for(int j =0; j<n_isotopes; j++){
      mutotal(i,j) = mutop(i,j)/mubtm(i,j);
    }
  }
  
  // Now do the same for sigma
  // 
  NumericMatrix sigsqtop(n, n_isotopes);
  NumericMatrix sigsqbtm(n, n_isotopes);
  
  for(int j=0; j<n_isotopes; j++){
    for(int i=0; i<n; i++){
      for(int k=0; k<n_sources; k++){
        sigsqtop(i,j) += pow(p(i,k) * concentrationmeans(k,j),2) * (pow(sourcesds(k,j),2) + pow(corrsds(k,j),2));
        sigsqbtm(i,j) += pow(p(i,k)*concentrationmeans(k,j),2);
      }
    }
  }
  
  
  
  NumericMatrix sigtotal(n,n_isotopes);
  
  // for(int i=0; i<(n); i++){
  //   for(int j=0; j<n_isotopes; j++){
  //     sigtotal(i,j) = pow((sigtop(i,j)/sigbtm(i,j)) + 1/theta(j + n_sources*n_covariates), 0.5);
  //   }
  // }
  
  for(int i=0; i<(n); i++){
    for(int j=0; j<n_isotopes; j++){
      
      // this is the sqrt of p^2 + q^2 etc + 1/tau
      sigtotal(i,j) = pow(sigsqtop(i,j)/sigsqbtm(i,j) + 1/theta(j + n_sources*n_covariates), 0.5);
    }
  }
  
  
  double hold = 0;
  
  // extract theta before below step
  // Then rewrite formula - drop n
  
  //double ly = y.rows();
  
  for(int i=0; i<n; i++){
    for(int j=0; j<n_isotopes; j++){
      hold += - log(sigtotal(i,j)) - 0.5  * log(2 * M_PI) -
        0.5 * pow((y(i,j)-mutotal(i,j)),2) * 1/(pow(sigtotal(i,j), 2));
    }
  }
  
  // for(int i=0; i<n; i++){
  //   for(int j=0; j<n_isotopes; j++){
  //     hold += - ly * log(sigtotal(i,j)) - 0.5 * ly * log(2 * M_PI) -
  //       0.5 * pow((y(i,j)-mutotal(i,j)),2) * 1/pow(sigtotal(i,j), 2);
  //   }
  // }
  
  
  
  double betanorm = 0;
  
  
  
  
  // This is what I've just added - was above
  // for(int i = 0; i<n_covariates; i++){
  //   for(int j=0; j<n_sources; j++){
  //     betanorm +=  - log(prior_sd(i,j)) - 0.5 * log(2 * M_PI) -
  //       0.5 * (pow((beta(i,j) - prior_means(i,j)), 2) * 1/(pow(prior_sd(i,j), 2)));
  //   }
  // }
  
  for(int i = 0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
      betanorm +=  -n_covariates * n_sources *  log(prior_sd(i,j)) - 0.5 * log(2 * M_PI) -
        0.5 * (pow((beta(i,j) - prior_means(i,j)), 2) * 1/(pow(prior_sd(i,j), 2)));
    }
  }
  
  
  double gammaprior = 0;
  
  for (int i=0; i<(n_isotopes); i++){
    gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +
      (c_0(i) - 1) * log(theta(i+n_sources*n_covariates))-
      d_0(i) * theta(i+n_sources*n_covariates);
    
  }
  
  // for (int i=0; i<(n_isotopes); i++){
  //   gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +
  //     (c_0(i) - 1) * (theta(i+n_sources*n_covariates))-
  //     d_0(i) * theta(i+n_sources*n_covariates);
  // 
  // }
  // 
  
  double totx = 0;
  totx = gammaprior + betanorm + hold;
  
  return (totx);
  
}


// This is basically the same as sim_theta but its using updated lambdas
// instead of set prior values

//[[Rcpp::export]]
double log_q_cpp(NumericVector theta, NumericVector lambda, 
                 int n_sources, int n_tracers, int S, int n_cov){
  
  NumericVector mean_beta(n_cov * n_sources);
  
  int mat_size = ((n_cov*n_sources) * (n_cov*n_sources+1))/2;
  
  
  for(int i=0; i<n_cov*n_sources; i++){
    mean_beta(i) = lambda(i);
  }
  
  
  
  NumericVector sig_beta(mat_size);
  
  
  for(int i = 0; i<mat_size; i++){
    sig_beta(i) = lambda(i + n_cov*n_sources);
    
  }
  
  
  
  NumericMatrix chol_prec(n_cov * n_sources, n_cov * n_sources);
  
  int count = 0;
  
  for(int i = 0; i<n_sources* n_cov; i++){
    for(int j = 0; j<n_sources * n_cov; j++){
      if (i <= j){
        count +=1;
        chol_prec((i),(j)) = sig_beta(count-1);
        
        
      }
      else{
        chol_prec(i,j) = 0;
      }
      
    }
  }
  
  NumericMatrix prec(n_sources*n_cov, n_sources*n_cov);
  prec = crossprod(chol_prec);
  
  
  NumericMatrix betaminusmean(1, n_sources * n_cov);
  
  
  for(int i = 0; i<n_sources*n_cov; i++){
    
    betaminusmean(0,i) = theta(i) - lambda(i);
  }
  
  NumericMatrix Z(1, n_sources*n_cov);
  
  NumericMatrix tcholprec(n_sources*n_cov, n_sources*n_cov);
  
  for(int i=0; i<n_sources*n_cov; i++){
    for(int j=0; j<n_sources*n_cov; j++){
      tcholprec(i,j) = chol_prec(j,i);
    }
  }
  
  
  Z = matmult(betaminusmean, tcholprec);
  
  NumericMatrix tZ(n_sources*n_cov, 1);
  
  
  for(int i = 0; i<n_sources*n_cov; i++){
    tZ(i,0) = Z(0,i);
  }
  
  
  NumericMatrix ZtZmat(1,1);
  
  ZtZmat = matmult(Z, tZ);
  
  double ZtZ = ZtZmat(0,0);
  
  double prod_sig;
  
  NumericVector sig(n_cov*n_sources);
  
  
  prod_sig = log(proddiag(tcholprec));
  
  
  double thetanorm = 0;
  
  thetanorm = -  (n_sources*n_cov/2) * log(2 * M_PI) + prod_sig - 0.5 * ZtZ;
  
  
  
  
  
  double gamman = 0;
  for (int i=0; i <(n_tracers); i++){
    gamman += lambda(mat_size + n_cov *n_sources +i) *
      log(lambda(mat_size + n_cov *n_sources +i + n_tracers))  -
      log(tgamma(lambda(mat_size + n_cov *n_sources +i)))  +
      (lambda(mat_size + n_cov *n_sources +i) - 1) * log(theta(i+n_sources*n_cov)) -
      lambda(mat_size + n_cov *n_sources +i + n_tracers) * theta((i+n_sources*n_cov));
  }
  
  // for (int i=0; i <(n_tracers); i++){
  //   gamman += lambda(n_covariates * mat_size + n_covariates *n_sources +i) *
  //     log(lambda(n_covariates * mat_size + n_covariates *n_sources +i + n_tracers))  -
  //     log(tgamma(lambda(n_covariates * mat_size + n_covariates *n_sources +i)))  +
  //     (lambda(n_covariates * mat_size + n_covariates *n_sources +i) - 1) * (theta(i+n_sources*n_covariates)) -
  //     lambda(n_covariates * mat_size + n_covariates *n_sources +i + n_tracers) * theta((i+n_sources*n_covariates));
  // }
  
  
  
  
  //double x = sum_p+ gamman;
  double x = thetanorm +gamman;
  
  return (x);
  
}

// TRY ADDING AUTO-DIFFERRENTIATION HERE RATHER THAN FROM FIRST PRINCIPLES
// This is just differentiating from first principles

// [[Rcpp::export]]
NumericVector delta_lqltcpp(NumericVector lambda, NumericVector theta, 
                            double eps, int n_sources, 
                            int n_tracers, int n_covariates,
                            int S) {
  // eps = 0.001;
  double k = lambda.length();
  NumericVector ans(k);
  NumericVector d(k);
  NumericVector lambdaplusd(k);
  NumericVector lambdaminusd(k);
  
  
  
  for(int i = 0; i<k; i++){
    
    for (int j = 0; j<k; j++){
      d(j) = 0;
    }
    d(i) = eps;
    
    
    for (int j = 0; j<k; j++){
      lambdaplusd(j) = lambda(j) + d(j);
      lambdaminusd(j) = lambda(j) - d(j);
    }
    ans(i) = (log_q_cpp(theta, lambdaplusd, n_sources, n_tracers, S, n_covariates) -  
      log_q_cpp(theta, lambdaminusd, n_sources, n_tracers, S, n_covariates))/(2 * eps);
  }
  return  ans;
}

// getting the difference of hcpp and log_q_cpp
// [[Rcpp::export]]
double h_lambdacpp(int n_sources, int n_isotopes,
                   NumericVector beta_prior,
                   int n_covariates,
                   int S,
                   NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                   NumericMatrix correctionmeans,
                   NumericMatrix corrsds, NumericMatrix sourcesds,
                   NumericVector theta, NumericMatrix y,
                   NumericVector lambda,
                   NumericMatrix x_scaled,
                   NumericVector c_0) {
  
  return hcpp(n_sources, n_isotopes, n_covariates, beta_prior, x_scaled, concentrationmeans, sourcemeans, correctionmeans,
              corrsds, sourcesds, theta, y, c_0) - log_q_cpp(theta, lambda, n_sources, n_isotopes, S, n_covariates);
}

//calculating covariance of 2 matrices
// [[Rcpp::export]]
NumericMatrix cov_mat_cpp(NumericMatrix x, NumericMatrix y) {
  int xcol = x.ncol();
  int ycol = y.ncol();
  int xrow = x.nrow();
  int yrow = y.nrow();
  
  NumericVector meanx(xcol);
  NumericVector meany(ycol); 
  NumericMatrix covmat(xcol, ycol);
  
  
  for(int i = 0; i<xcol; i++){
    meanx(i) = mean(x(_,i));
  }
  for(int i = 0; i<ycol; i++){
    meany(i) = mean(y(_,i));
  }
  
  NumericMatrix xminusmean(xrow, xcol);
  NumericMatrix yminusmean(yrow, ycol);
  
  for(int j = 0; j<xcol; j++){
    for(int i=0; i<xrow; i++){
      xminusmean(i,j) = x(i,j) - meanx(j);
    }
  }
  
  for(int j = 0; j<ycol; j++){
    for(int i =0; i<yrow; i++){
      yminusmean(i,j) = y(i,j) - meany(j);
    }
  }
  
  NumericMatrix sumxy(xcol, ycol);
  
  // NumericVector xcol(x.ncol());
  // NumericVector ycol(y.ncol());
  
  for(int i = 0; i<xcol; i++){
    for(int j=0; j<ycol; j++){
      for(int n =0; n<xrow; n++){
        
        sumxy(i,j) += xminusmean(n,i) * yminusmean(n,j);
      }
    }}
  
  
  for(int i=0; i<xcol; i++){
    for(int j = 0; j<ycol; j++){
      covmat(i,j) = sumxy(i,j)/(xrow-1);
    }
  }
  
  return covmat;
}


//Nabla LB is the mean of delta_lqlt element-wise multiplied by h_lambda
// [[Rcpp::export]]
NumericVector nabla_LB_cpp(NumericVector lambda, NumericMatrix theta, 
                           int n_sources, int n_tracers, NumericVector beta_prior,
                           int S, int n_covariates,
                           NumericMatrix x_scaled,
                           NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                           NumericMatrix correctionmeans,
                           NumericMatrix corrsds, NumericMatrix sourcesds, NumericMatrix y,
                           NumericVector c,
                           NumericVector c_0){
  
  int thetanrow = theta.nrow();
  int lambdalength = lambda.length();
  
  NumericMatrix big_c(thetanrow, c.length());
  
  //working
  NumericMatrix big_delta_lqlt(thetanrow, lambdalength); 
  NumericMatrix big_h_lambda_rep(lambdalength, thetanrow);
  NumericMatrix big_h_lambda_rep_transpose(thetanrow, lambdalength);
  
  NumericVector big_h_lambda(thetanrow);
  NumericVector big_h_lambda_transpose(thetanrow);
  
  
  for(int i = 0; i <thetanrow; i++){
    big_delta_lqlt(i,_) = delta_lqltcpp(lambda, theta(i,_), 0.001, n_sources, n_tracers,
                   n_covariates, S);
  }
  
  for(int i =0; i<thetanrow; i++){
    big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers, beta_prior,
                 n_covariates, S,
                 concentrationmeans, sourcemeans,
                 correctionmeans,
                 corrsds,sourcesds, theta(i,_), y,
                 lambda, x_scaled, c_0);
  }
  
  
  
  for(int i =0; i<lambdalength; i++){
    big_h_lambda_rep(i,_) = big_h_lambda;
  }
  
  for(int i=0; i<lambdalength; i++){
    for (int j=0; j < theta.nrow(); j++){
      big_h_lambda_rep_transpose(j,i) = big_h_lambda_rep(i,j);
    }}
  
  
  for(int i =0; i<thetanrow; i++){
    big_c(i,_) = c;
  }
  
  
  
  
  
  NumericMatrix big_h_minus_c(thetanrow, lambdalength);
  //NumericMatrix big_h_minus_c_t(lambda.length(), theta.nrow());
  
  for (int i = 0; i<thetanrow; i++){
    for(int j = 0; j<lambdalength; j++){
      big_h_minus_c(i,j) = big_h_lambda_rep_transpose(i,j) - big_c(i,j);
    }
  }
  
  //big_h_minus_c_t = transpose(big_h_minus_c);
  
  NumericMatrix ansmat(big_delta_lqlt.nrow(), big_h_minus_c.ncol());
  
  for (int i = 0; i < big_delta_lqlt.nrow(); i++) 
  {
    for (int j = 0; j < big_delta_lqlt.ncol(); j++) {
      
      
      ansmat(i,j) = big_delta_lqlt(i,j) * big_h_minus_c(i,j);
      
      
    }
  }
  
  NumericVector ans(ansmat.ncol());
  for(int i = 0; i<ansmat.ncol(); i++){
    
    ans(i) = mean(ansmat(_,i));
    
  }
  
  return ans;
}


// calculate control variate (big formula)
// [[Rcpp::export]]
NumericVector control_var_cpp(NumericVector lambda, 
                              NumericMatrix theta, 
                              int n_sources, int n_tracers,
                              NumericVector beta_prior,
                              int n_covariates,
                              NumericMatrix x_scaled,
                              NumericMatrix concentrationmeans, 
                              NumericMatrix sourcemeans,
                              NumericMatrix correctionmeans,
                              NumericMatrix corrsds, 
                              NumericMatrix sourcesds, 
                              NumericMatrix y,
                              NumericVector c_0){
  
  int S = theta.nrow();
  int lambdallength = lambda.length();
  NumericMatrix big_delta_lqlt(S, lambdallength); 
  NumericMatrix big_h_lambda_rep(lambdallength, S);
  NumericMatrix big_h_lambda_rep_transpose(S, lambdallength);
  NumericVector big_h_lambda(S);
  NumericVector big_h_lambda_transpose(S);
  
  for(int i = 0; i <S; i++){
    big_delta_lqlt(i,_) = delta_lqltcpp(lambda, theta(i,_), 0.001, n_sources, n_tracers,
                   n_covariates, S);
  }
  
  for(int i =0; i<S; i++){
    big_h_lambda(i) = h_lambdacpp(n_sources, n_tracers, beta_prior,
                 n_covariates, S,
                 concentrationmeans, sourcemeans,
                 correctionmeans,
                 corrsds,sourcesds, theta(i,_), y,
                 lambda, x_scaled, c_0);
  }
  
  
  for(int i =0; i<lambdallength; i++){
    big_h_lambda_rep(i,_) = big_h_lambda;
  }
  
  for(int i=0; i<lambdallength; i++){
    for (int j=0; j < theta.nrow(); j++){
      big_h_lambda_rep_transpose(j,i) = big_h_lambda_rep(i,j);
    }}
  
  NumericMatrix big_nabla(S, lambdallength);
  
  for (int i = 0; i < S; i++)
  {
    for (int j = 0; j < lambdallength; j++) {
      
      
      big_nabla(i,j) = big_delta_lqlt(i,j) * big_h_lambda_rep_transpose(i,j);
      
      
    }
  }
  
  NumericVector var_big_delta_lqlt(lambdallength);
  
  for(int i = 0; i<lambdallength; i++){
    var_big_delta_lqlt(i) = var(big_delta_lqlt(_,i));
  }
  
  NumericMatrix covmat(lambdallength, lambdallength);
  
  covmat = cov_mat_cpp(big_nabla, big_delta_lqlt);
  
  NumericVector diag(lambdallength);
  for(int i =0; i<lambdallength; i++){
    for(int j =0; j<lambdallength; j++){
      if(i == j){
        diag(i) = covmat(i,j);
      }
    }}
  
  NumericVector ans(lambdallength);
  for(int i =0; i<lambdallength; i++){
    ans(i) = diag(i)/var_big_delta_lqlt(i);
  }
  
  return ans;
}

// estimate of LB 
// [[Rcpp::export]]
double LB_lambda_cpp(NumericMatrix theta, NumericVector lambda, 
                     NumericVector p, int n_sources, int n_isotopes, 
                     NumericVector beta_prior,
                     int n_covariates,
                     NumericMatrix x_scaled,
                     NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
                     NumericMatrix correctionmeans,
                     NumericMatrix corrsds, 
                     NumericMatrix sourcesds, 
                     NumericMatrix y,
                     NumericVector c_0){
  int S = theta.nrow();
  
  NumericVector hlambdaapply(S);
  
  for(int i = 0; i <S; i++){
    hlambdaapply(i) = h_lambdacpp(n_sources, n_isotopes, beta_prior, 
                 n_covariates, S,
                 concentrationmeans, sourcemeans,
                 correctionmeans, corrsds, sourcesds, 
                 theta(i,_), y, lambda,
                 x_scaled, c_0);
  }
  
  double ans = mean(hlambdaapply);
  
  return ans;
  
  
}


// Actually putting it all together and running it
// [[Rcpp::export]]
NumericVector run_VB_cpp(NumericVector lambdastart,
                         int n_sources,
                         int n_tracers,
                         int n_covariates,
                         int n,
                         NumericVector beta_prior,
                         NumericMatrix concentrationmeans,
                         NumericMatrix sourcemeans,
                         NumericMatrix correctionmeans,
                         NumericMatrix corrsds,
                         NumericMatrix sourcesds,
                         NumericMatrix y,
                         NumericMatrix x_scaled,
                         int S,
                         int P,
                         double beta_1,
                         double beta_2,
                         int tau,
                         double eps_0,
                         int t_W,
                         NumericVector c_prior,
                         bool solo
){
  
  
  
  int lsl = lambdastart.length();
  
  NumericMatrix theta(S, (n_sources + n_tracers));
  
  theta = sim_thetacpp(S, lambdastart, n_sources, n_tracers, n_covariates, solo);
  // print(theta);
  NumericVector c(lsl);
  
  c = control_var_cpp(lambdastart, theta, n_sources, n_tracers, beta_prior, 
                      n_covariates, x_scaled,
                      concentrationmeans,
                      sourcemeans, correctionmeans,
                      corrsds, sourcesds, y, c_prior);
  
  NumericVector g_0(lsl);
  
  NumericVector c_0(lsl);
  for(int i=0; i<lsl; i++){
    c_0(i) = 0;
  }
  
  g_0 = nabla_LB_cpp(lambdastart, theta, 
                     n_sources, n_tracers, 
                     beta_prior,
                     S, n_covariates, x_scaled,
                     concentrationmeans, sourcemeans,
                     correctionmeans, corrsds, 
                     sourcesds, y, c_0, c_prior);
  
  NumericVector nu_0(lsl);
  
  for(int i = 0; i<lsl; i++){
    nu_0(i) = pow(g_0(i),2);
  }
  
  
  NumericVector g_bar(lsl);
  g_bar = g_0;
  
  NumericVector nu_bar(lsl);
  nu_bar = nu_0;
  
  NumericVector g_t(lsl);
  NumericVector nu_t(lsl);
  
  double patience = 0;
  bool stop = FALSE;
  double max_LB_bar = R_NegInf;
  double alpha_t = 0;
  double t = 0;
  NumericVector LB(t_W+1);
  
  for(int i=0; i<t_W; i++){
    LB(i) = NA_REAL;
  }
  
  NumericVector lambda(lsl);
  
  for(int i = 0; i<lsl; i++){
    lambda(i) = lambdastart(i);
  }
  
  while(stop == FALSE){
    
    theta = sim_thetacpp(S, lambda, n_sources, n_tracers, n_covariates, solo);
    
    
    g_t = nabla_LB_cpp(lambda, theta, n_sources, n_tracers, beta_prior,
                       S, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds, sourcesds,
                       y, c, c_prior);
    
    
    c = control_var_cpp(lambda, theta,n_sources,n_tracers, beta_prior,
                        n_covariates, x_scaled,
                        concentrationmeans, sourcemeans,
                        correctionmeans,
                        corrsds,sourcesds, y, c_prior);
    
    
    for(int i=0; i<lsl; i++){
      nu_t(i) = pow(g_t(i),2);
    }
    
    
    for(int i=0; i<lsl; i++){
      g_bar(i) = (beta_1 * g_bar(i)) + ((1-beta_1) * g_t(i));
      nu_bar(i) = (beta_2 * nu_bar(i)) + ((1-beta_2) * nu_t(i));
    }
    
    NumericVector alpha_min(2);
    
    for(int i=0; i<2; i++){
      alpha_min(0) = eps_0;
      alpha_min(1) = eps_0 * (tau/(t+1));
    }
    
    alpha_t = Rcpp::min(alpha_min);
    
    
    //# Update lambda
    for(int i = 0; i<lsl; i++){
      //lambda(i) = lambda(i) + alpha_t * 1/pow(nu_bar(i), 0.5);
      lambda(i) = lambda(i) + alpha_t * (g_bar(i)/(pow(nu_bar(i), 0.5)));
      
      //Adding if statement to stop lambda dropping below 0 - not sure if this
      // is 1. allowed 2. will mess up results 3. will result in it just converging
      // on 0.001. Had it smaller but it broke it. Actually it just seems to break it
      //full stop. Or - it means the issue isn't with lambda going below zero?
      //if(lambda(lsl-1) < 0){lambda(lsl-1) = 0.001;}
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
    
    //# Compute the moving average LB if out of warm-up
    if(t<=t_W){
      
      // # Compute a new lower bound estimate
      NumericVector p = hfn(theta, n_sources, n,  n_covariates, x_scaled);
      
      LB(t) = LB_lambda_cpp(theta,lambda, p, n_sources, n_tracers,
         beta_prior,
         n_covariates, x_scaled,
         concentrationmeans, sourcemeans,
         correctionmeans,
         corrsds,sourcesds, y, c_prior);
    }
    else{
      for (int i = 0; i<(t_W-1); i++){
        LB(i) = LB(i+1);
      }
      
      NumericVector p = hfn(theta, n_sources, n, n_covariates, x_scaled);
      
      LB(t_W) = LB_lambda_cpp(theta, lambda, p, n_sources, n_tracers,
         beta_prior,
         n_covariates, x_scaled,
         concentrationmeans, sourcemeans,
         correctionmeans,
         corrsds,sourcesds, y, c_prior);
      
      double LB_bar = mean(LB);
      
      NumericVector maxbar(2);
      
      for(int i=0; i<2; i++){
        maxbar(0) = max_LB_bar;
        maxbar(1) = LB_bar;
      }
      
      max_LB_bar = Rcpp::max(maxbar);
      
      
      if(LB_bar>= max_LB_bar){
        patience = 0;
      } else{
        patience = patience +1;
      }
      
    }
    
    //////////////////////
    if(patience>P){
      stop = TRUE;
    }
    
    if(t>500){
      stop = TRUE;
    }
    t = t + 1;
  }
  
  
  return lambda;
  
}
