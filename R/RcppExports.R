# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

digamma_wrapper <- function(x) {
    .Call('_cosimmr_digamma_wrapper', PACKAGE = 'cosimmr', x)
}

rMVNormCpp <- function(n, Mean, Var) {
    .Call('_cosimmr_rMVNormCpp', PACKAGE = 'cosimmr', n, Mean, Var)
}

sim_thetacpp <- function(S, lambda, n_sources, n_tracers, n_cov, solo) {
    .Call('_cosimmr_sim_thetacpp', PACKAGE = 'cosimmr', S, lambda, n_sources, n_tracers, n_cov, solo)
}

hfn <- function(theta, n_sources, n, n_cov, x_scaled) {
    .Call('_cosimmr_hfn', PACKAGE = 'cosimmr', theta, n_sources, n, n_cov, x_scaled)
}

hcpp <- function(n_sources, n_isotopes, n_covariates, d_prior, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, theta, y, c_prior, sd_prior) {
    .Call('_cosimmr_hcpp', PACKAGE = 'cosimmr', n_sources, n_isotopes, n_covariates, d_prior, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, theta, y, c_prior, sd_prior)
}

log_q_cpp <- function(theta, lambda, n_sources, n_tracers, S, n_cov) {
    .Call('_cosimmr_log_q_cpp', PACKAGE = 'cosimmr', theta, lambda, n_sources, n_tracers, S, n_cov)
}

delta_lqltcppauto <- function(lambda, theta, n_sources, n_tracers, n_covariates, S) {
    .Call('_cosimmr_delta_lqltcppauto', PACKAGE = 'cosimmr', lambda, theta, n_sources, n_tracers, n_covariates, S)
}

h_lambdacpp <- function(n_sources, n_isotopes, beta_prior, n_covariates, S, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, theta, y, lambda, x_scaled, c_0, sd_prior) {
    .Call('_cosimmr_h_lambdacpp', PACKAGE = 'cosimmr', n_sources, n_isotopes, beta_prior, n_covariates, S, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, theta, y, lambda, x_scaled, c_0, sd_prior)
}

cov_mat_cpp <- function(x_arma, y_arma) {
    .Call('_cosimmr_cov_mat_cpp', PACKAGE = 'cosimmr', x_arma, y_arma)
}

nabla_LB_cpp <- function(lambda, theta, n_sources, n_tracers, beta_prior, S, n_covariates, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, c, c_0, sd_prior) {
    .Call('_cosimmr_nabla_LB_cpp', PACKAGE = 'cosimmr', lambda, theta, n_sources, n_tracers, beta_prior, S, n_covariates, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, c, c_0, sd_prior)
}

control_var_cpp <- function(lambda, theta, n_sources, n_tracers, beta_prior, n_covariates, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, c_0, sd_prior) {
    .Call('_cosimmr_control_var_cpp', PACKAGE = 'cosimmr', lambda, theta, n_sources, n_tracers, beta_prior, n_covariates, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, c_0, sd_prior)
}

LB_lambda_cpp <- function(theta, lambda, n_sources, n_isotopes, beta_prior, n_covariates, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, c_0, sd_prior) {
    .Call('_cosimmr_LB_lambda_cpp', PACKAGE = 'cosimmr', theta, lambda, n_sources, n_isotopes, beta_prior, n_covariates, x_scaled, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, c_0, sd_prior)
}

run_VB_cpp <- function(lambdastart, n_sources, n_tracers, n_covariates, n, beta_prior, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, x_scaled, S, P, beta_1, beta_2, tau, eps_0, t_W, c_prior, solo, sd_prior) {
    .Call('_cosimmr_run_VB_cpp', PACKAGE = 'cosimmr', lambdastart, n_sources, n_tracers, n_covariates, n, beta_prior, concentrationmeans, sourcemeans, correctionmeans, corrsds, sourcesds, y, x_scaled, S, P, beta_1, beta_2, tau, eps_0, t_W, c_prior, solo, sd_prior)
}

