#' Run a \code{cosimmr_input} object through STAN 
#'
#' This is the main function of cosimmr. It takes a \code{cosimmr_input} object
#' created via \code{\link{cosimmr_load}}, runs it in fixed form
#' Variational Bayes  via STAN to determine the dietary proportions, and then
#' outputs a \code{cosimmr_output} object for further analysis and plotting
#' via \code{\link{plot.cosimmr_output}}.
#'
#'
#'@param type What type of model to run using STAN. Options are 'STAN_VB
#'@param error_type Whether to use 'processxresidual' error term or 
#''process+residual' term. Defaults to 'processxresidual'
#'@param prior_control A list of values including arguments named \code{sigma_shape} 
#'(prior values for sigma shape), \code{sigma_rate} (prior values for sigma rate)
#'@param n_samples Number of samples to output. Defaults to 3600.
#'
#'@return an object of class \code{cosimmr_output} with two named top-level 
#'components: \item{input }{The \code{cosimmr_input} object given to the
#' \code{cosimmr_ffvb} function} \item{output }{A set of outputs produced by
#' the FFVB function. These can be analysed using the
#' \code{\link{summary.cosimmr_output}} and \code{\link{plot.cosimmr_output}}
#' functions.}
#'
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
#'
#' @seealso \code{\link{cosimmr_load}} for creating objects suitable for this
#' function, \code{\link{plot.cosimmr_input}} for creating isospace plots,
#' \code{\link{summary.cosimmr_output}} for summarising output, and
#' \code{\link{plot.cosimmr_output}} for plotting output.
#'
#' @references Andrew C. Parnell, Donald L. Phillips, Stuart Bearhop, Brice X.
#' Semmens, Eric J. Ward, Jonathan W. Moore, Andrew L. Jackson, Jonathan Grey,
#' David J. Kelly, and Richard Inger. Bayesian stable isotope mixing models.
#' Environmetrics, 24(6):387–399, 2013.
#'
#' Andrew C Parnell, Richard Inger, Stuart Bearhop, and Andrew L Jackson.
#' Source partitioning using stable isotopes: coping with too much variation.
#' PLoS ONE, 5(3):5, 2010.
#'
#'
#' @examples
#' \donttest{
#' ## See the package vignette for a detailed run through of these examples
#'
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' data(geese_data_day1)
#' x =  c(1,2,3,2,1,3,2,1,2)
#' cosimmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ x,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means
#'   )
#' )
#'
#' # Plot
#' plot(cosimmr_1)
#'
#' # Print
#' cosimmr_1
#'
#' # FFVB run
#' cosimmr_1_out <- cosimmr_stan(cosimmr_1)
#'
#' # Print it
#' print(cosimmr_1_out)
#'
#' # Summary
#' summary(cosimmr_1_out, type = "correlations")
#' summary(cosimmr_1_out, type = "statistics")
#' ans <- summary(cosimmr_1_out, type = c("quantiles", "statistics"))
#'
#' # Plot
#' plot(cosimmr_1_out, type = "beta_boxplot")
#' plot(cosimmr_1_out, type = "beta_histogram")
#'
#'}
#' @export cosimmr_ffvb
cosimmr_stan <- function(cosimmr_in,
                         type = "STAN_VB",
                         error_type = "processxresidual",
                         prior_control = list(
                          sigma_shape = c(rep(1, cosimmr_in$n_tracers)),
                          sigma_rate = c(rep(1, cosimmr_in$n_tracers))
                         ),
                         n_samples = 3600
                         
                         ){
  
  #Core detection - potentially have this as an option people can turn on and off?
  options(mc.cores = parallel::detectCores())

  #Need loops for error type
  if(error_type == "process+residual"){


  #Need to add solo check here
  if (nrow(cosimmr_in$mixtures) == 1) {
    message("Only 1 mixture value, performing a simmr solo run...\n")
    solo <- 0
  } else {
    solo <- 1
  }


  stan_dat = list(
    J = cosimmr_in$n_tracers,
    N = cosimmr_in$n_obs,
    K = cosimmr_in$n_sources,
    L = cosimmr_in$n_covariates,
    y = cosimmr_in$mixtures,
    q = cosimmr_in$concentration_means,
    s_mean = cosimmr_in$source_means,
    c_mean = cosimmr_in$correction_means,
    s_sd = cosimmr_in$source_sds,
    c_sd = cosimmr_in$correction_sds,
    sigma_shape = prior_control$sigma_shape,
    sigma_rate = prior_control$sigma_rate,
    X = cosimmr_in$x_scaled,
    not_solo = solo
  )

  if(type == "STAN_VB"){
    #I think this is calling the model - need to check that it calls what I actually need
    #Don't need this - model is compiled already afaik
   # model = stan_model(stanmodels$STAN_VB)

  # Fit using VB
  # Get good starting value sby optimizing first

  fit_opt <- optimizing(stanmodels$STAN_VB,
                        data = stan_dat)
  
  which_beta <- grep('beta', names(fit_opt$par))
  
  which_sigma <- grep('sigma', names(fit_opt$par))
  
  beta_start_opt <- structure(fit_opt$par[which_beta],
                              dim  = c(stan_dat$K, stan_dat$L))
  sigma_raw_start_opt <- fit_opt$par[which_sigma][1:2]
  
  
  fit_vb <- vb(
    stanmodels$STAN_VB, data = stan_dat,
    algorithm = 'fullrank',
    pars = c('beta', 'sigma'),
    init = list('beta' = beta_start_opt,
                'log_sigma_raw' = sigma_raw_start_opt),
    tol_rel_obj = 0.00001, #convergence tolerance on the relative norm of the objective
    output_samples = n_samples
     )
  
  #Want to simulate p here and return that?
  #And return sigma in a sensible way
  extracted_samples = rstan::extract(fit_vb)

  #Want to extract all the betas in a sensible way first I think
  beta_ans = extracted_samples$beta # This is n_samples * K * n_covariates
  sigma_ans = extracted_samples$sigma
  
  
  #Want to simulate p
  

  mylist <- list(
fit_vb = fit_vb,
beta = beta_ans,
sigma = sigma_ans
  )
  
  
  output_all <- list(input = cosimmr_in, output = mylist)
  
  class(output_all) <- c("cosimmr_output", "stan_vb")
  
  return(output_all)

  


  } else if(type == "STAN_MCMC"){
    #I think this is calling the model - need to check that it calls what I actually need
    model = stan_model(stanmodels$STAN_VB)
    fit_mcmc <- sampling(
      model,
      data = stan_dat,
      seed = 1
    )
  }
  } else if(error_type == "processxresidual"){


  }
  
}



# # Fit using VB
# # Get good starting value sby optimizing first
# fit_opt <- optimizing(model, 
#                       data = stan_dat)
# which_beta <- grep('beta', names(fit_opt$par))
# which_sigma <- grep('sigma', names(fit_opt$par))
# beta_start_opt <- structure(fit_opt$par[which_beta],
#                             dim  = c(K, L))
# sigma_raw_start_opt <- fit_opt$par[which_sigma][1:2]
# # beta_start <- structure(c(0.963, -1.119, -0.403, 0.457, -0.099, 0.364, -1.085,
# #                           0.362), dim = c(4L, 2L))
# # log_sigma_raw_start <- c(1.5, 0.02)
# fit_vb <- vb(
#   model, data = stan_dat, 
#   algorithm = 'fullrank',
#   pars = c('beta', 'sigma'),
#   init = list('beta' = beta_start_opt, 
#               'log_sigma_raw' = sigma_raw_start_opt), 
#   tol_rel_obj = 0.00001,
#   seed = 123
# )
# fit_vb
# plot(fit_vb, par = 'beta')
# plot(fit_vb, par = 'sigma')
# stop()
# 
# # Try just finding best values
# 
# # Fit using MCMC
# fit_mcmc <- sampling(
#   model, 
#   data = stan_dat, 
#   seed = 1
# )
# fit_mcmc
# plot(fit_mcmc, par = 'beta')
# plot(fit_mcmc, par = 'sigma')
# 
# 
# # Compare with simmr to check it looks the same - it does)
# simmr_1 <- with(
#   geese_data,
#   simmr_load(
#     mixtures = mixtures,
#     source_names = source_names,
#     source_means = source_means,
#     source_sds = source_sds,
#     correction_means = correction_means,
#     correction_sds = correction_sds,
#     concentration_means = concentration_means
#   )
# )
# 
# # MCMC run
# simmr_1_out <- simmr_mcmc(simmr_1)
# 
# # Print it
# summary(simmr_1_out)
