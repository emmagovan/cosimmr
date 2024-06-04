#' Run a \code{cosimmr_input} object through the Fixed Form Variational
#' Bayes(FFVB) function
#'
#' This is the main function of cosimmr. It takes a \code{cosimmr_input} object
#' created via \code{\link{cosimmr_load}}, runs it in fixed form
#' Variational Bayes to determine the dietary proportions, and then
#' outputs a \code{cosimmr_output} object for further analysis and plotting
#' via \code{\link{plot.cosimmr_output}}.
#'

#' @param cosimmr_in An object created via the function \code{\link{cosimmr_load}}
#' @param prior_control A list of values including arguments named \code{mu_0}
#' (prior for mu), and \code{sigma_0} (prior for sigma).
#' @param ffvb_control A list of values including arguments named \code{n_output}
#' (number of rows in theta output), \code{S} (number of samples taken at each
#' iteration of the algorithm), \code{P} (patience parameter), \code{beta_1}
#' and \code{beta_2} (adaptive learning weights), \code{tau} (threshold for
#' exploring learning space), \code{eps_0} (fixed learning rate),
#' \code{t_W} (rolling window size)
#' @param error_type Option to choose "processxresidual" or "process+residual"
#' @return An object of class \code{cosimmr_output} with two named top-level
#' components: \item{input }{The \code{cosimmr_input} object given to the
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
#' @importFrom R2jags jags
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
#' cosimmr_1_out <- cosimmr_ffvb(cosimmr_1)
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
#' plot(cosimmr_1_out, type = "beta_boxplot", cov_name = "x")
#' plot(cosimmr_1_out, type = "beta_histogram", cov_name = "x")
#'
#'}
#' @export cosimmr_ffvb
cosimmr_ffvb <- function(cosimmr_in,
                         prior_control = list(
                           mu_0 = rep(0, (cosimmr_in$n_sources * cosimmr_in$n_covariates)), #this is currently both the prior value AND the starting lambda value - should change?
                           mu_log_sig_sq_0 = rep(0, cosimmr_in$n_tracers), #starting value for sigma mean
                           mu_log_sig_omi_0 = rep(1, cosimmr_in$n_tracers), #starting value for omicron mean
                           sigma_0 = 1, #repeated starting value for all sig and omicron
                           tau_shape = rep(1, cosimmr_in$n_tracers),
                           tau_rate = rep(1, cosimmr_in$n_tracers),
                           uniform_a = 0,
                           uniform_b = 20
                         ),
                         ffvb_control = list(
                           n_output = 3600, #3600
                           S = 500, #500
                           P = 50, #100
                           beta_1 = 0.75, #0.9
                           beta_2 = 0.75, #0.9
                           tau = 500, #50
                           eps_0 = 0.0011, #0.001
                           t_W =500 #1000
                         ),
                         error_type = "processxresidual") {
  # Throw a warning if less than 4 observations in a group - 1 is ok as it wil do a solo run
  #  if (min(table(cosimmr_in$group)) > 1 & min(table(cosimmr_in$group)) < 4) warning("At least 1 group has less than 4 observations - either put each observation in an individual group or use informative prior information")
  
  ## Just throw in a if else to make this two separate functions??
  if(error_type == "processxresidual"){
    output <- list #vector("list", length = cosimmr_in$n_groups)
    # names(output) <- levels(cosimmr_in$group)
    K <- cosimmr_in$n_sources
    n_tracers <- cosimmr_in$n_tracers
    n_output <- ffvb_control$n_output
    mu_a <- prior_control$mu_0
    sigma_a <- prior_control$sigma_0
    n_covariates<- cosimmr_in$n_covariates
    x_scaled = cosimmr_in$x_scaled
    uniform_a = prior_control$uniform_a
    uniform_b = prior_control$uniform_b
    
  
    
    #Regular
    beta_lambda <-c(mu_a, prior_control$mu_log_sig_sq_0,   prior_control$mu_log_sig_omi_0, rep(sigma_a, (((K*n_covariates +n_tracers*2) * ((K*n_covariates +n_tracers*2) +1))/2)))
    

    lambda <- c(beta_lambda)
    # ,
    #             prior_control$tau_shape, 
    #             prior_control$tau_rate)
    
    
    ll = length(lambda)
    
    lambdares <- c(rep(NA, ll))
    
    
    #Theta is all the betas (so K * n_obs * n_cov + 2*n_tracers)
    
    # thetares <- matrix(rep(NA, ((K * n_covariates + n_tracers* 2 * cosimmr_in$n_obs) * n_output)),
    #                    ncol = (K * n_covariates + n_tracers * 2 * cosimmr_in$n_obs),
    #                    nrow = n_output
    # )
    
    mylist <- list#vector("list", length = cosimmr_in$n_groups)
    
    # names(mylist) <- levels(cosimmr_in$group)
    
    p_fun <- function(x) exp(x) / sum(exp(x))

    
    # Determine if a single observation or not
    if (nrow(cosimmr_in$mixtures) == 1) {
      message("Only 1 mixture value, performing a simmr solo run...\n")
      solo <- TRUE
      beta_prior = c(rep(100, n_tracers))
    } else {
      solo <- FALSE
      beta_prior = prior_control$tau_rate
    }
    
    n_tracers <- cosimmr_in$n_tracers
    n_sources <- cosimmr_in$n_sources
    s_names <- cosimmr_in$source_names
    K <- cosimmr_in$n_sources
    source_means <- cosimmr_in$source_means
    source_sds <- cosimmr_in$source_sds
    correction_means <- cosimmr_in$correction_means
    correction_sds <- cosimmr_in$correction_sds
    concentration_means <- cosimmr_in$concentration_means
    y <- cosimmr_in$mixtures
    n = nrow(y)
    
    lambdaprior <- c(lambda)
    
    

    lambdastart = lambdaprior
   
    l_l = length(lambdastart)


    mu_beta_prior = matrix(prior_control$mu_0, nrow = cosimmr_in$n_covariates)
    sigma_beta_prior = matrix(c(rep(prior_control$sigma_0, cosimmr_in$n_covariates * cosimmr_in$n_sources)), nrow = cosimmr_in$n_covariates)
    
    
    lambdares <- run_VB_cpp_pxr(
      lambdastart, K, n_tracers, n_covariates, n, beta_prior, 
      concentration_means,
      source_means, correction_means, correction_sds,
      source_sds, y, (x_scaled), ffvb_control$S,
      ffvb_control$P, ffvb_control$beta_1,
      ffvb_control$beta_2, ffvb_control$tau,
      ffvb_control$eps_0, ffvb_control$t_W,
      prior_control$tau_shape, solo, sigma_beta_prior, mu_beta_prior,
      uniform_a, uniform_b)
    
    
    iteration = (lambdares$iteration) 
    a = which.max(lambdares$mean_LB_bar[(ffvb_control$t_W +2):iteration]) # +2 here for same reason, start from 0 and then greater than
    use_this = lambdares$lambda_save[(ffvb_control$t_W +1)+a,]

    kappa2 = rMVNormCpp(K*n_covariates + n_tracers*2, c(rep(0, n_output)), diag(n_output))
    
    thetares <- sim_thetacpp_pxr(n_output, use_this, K, n_tracers, n_covariates, solo, kappa2)

    sigma <- sqrt(exp((thetares[,(K*n_covariates + 1):(K*n_covariates + n_tracers)])))
    
    omicron <- 1/(1+exp(-(thetares[,(K*n_covariates + n_tracers + 1):(K*n_covariates + n_tracers*2)])))
    
    p_sample = array(NA, dim = c(cosimmr_in$n_obs, n_output, K))
    p_mean_sample = matrix(NA, nrow = n_output, ncol = K)

    beta = thetares[,1:(n_covariates * K)]
    
    f <- array(NA, dim = c(cosimmr_in$n_obs, K, n_output)) 

    f_mean_sample = matrix(NA, nrow = K, ncol = n_output)
    
    if(cosimmr_in$intercept == TRUE){
      x_sample_mean = c(1, rep(0, (ncol(cosimmr_in$x_scaled) - 1)))
    } else if(cosimmr_in$intercept == FALSE){
      x_sample_mean = c(rep(0, (ncol(cosimmr_in$x_scaled))))
    }
    
    for(s in 1:n_output){
      f[,,s] = as.matrix(cosimmr_in$x_scaled) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE)
      f_mean_sample[,s] = matrix(x_sample_mean, nrow = 1) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE) 
    }
    
    for(j in 1:n_output){
      for (n_obs in 1:cosimmr_in$n_obs) {
        p_sample[n_obs,j, ] <- exp(f[n_obs,1:K, j]) / (sum((exp(f[n_obs,1:K, j]))))
      }
    }
    
    for(j in 1:n_output){
      p_mean_sample[j,] = exp(f_mean_sample[1:K,j]) /(sum(exp(f_mean_sample[1:K,j])))
    }
    
    
    mylist <- list(
      source_names = cosimmr_in$source_names,
      theta = thetares,
      groupnames = cosimmr_in$group,
      lambdares = lambdares,#$lambda,
      # mean_LB = lambdares$mean_LB,
      beta = beta,
      BUGSoutput = list(
        sims.list = list(
          p = p_sample,
          p_mean = p_mean_sample,
          sigma = sigma,
          omicron = omicron
        )),
      model = list(data = list(
        mu_f_mean = c(mu_a),
        sigma_f_sd = c(rep(sigma_a, K)),
        ffvb_control = ffvb_control,
        prior_control = prior_control
      ))
    )
    
    
    output_all <- list(input = cosimmr_in, output = mylist)
    
    class(output_all) <- c("cosimmr_output", "ffvb")
    
    return(output_all)
    
  }
  else if(error_type == "process+residual"){
  output <- list #vector("list", length = cosimmr_in$n_groups)
  # names(output) <- levels(cosimmr_in$group)
  K <- cosimmr_in$n_sources
  n_tracers <- cosimmr_in$n_tracers
  n_output <- ffvb_control$n_output
  mu_a <- prior_control$mu_0
  sigma_a <- prior_control$sigma_0
  n_covariates<- cosimmr_in$n_covariates
  x_scaled = cosimmr_in$x_scaled
  
  # sig_prior = matrix(rep(0, n_covariates*K*n_covariates*K), ncol = n_covariates*K, 
  #                    nrow = n_covariates*K)
  # count = 0
  # sig_prior_vars = rep(sigma_a, (((K*n_covariates) * (K*n_covariates +1))/2))
  # 
  # for(i in 1:(n_covariates * K)){
  #   for(j in 1:(n_covariates * K)){
  #   if(j<=i){
  #     count = count +1
  #     sig_prior[i,j] = sig_prior_vars[count]
  #   }
  #   }
  # }
  
  # sig_prior = matrix(rep(sigma_a, n_covariates*K*n_covariates*K), ncol = n_covariates*K, 
  #                    nrow = n_covariates*K)
  
  # sig_prior = matrix(rep(sigma_a, (n_covariates*K+n_tracers)*(n_covariates*K+n_tracers)), ncol = (n_covariates*K + n_tracers), 
  #                    nrow = (n_covariates*K + n_tracers))
  
  #c_0 <- c(1,1)#prior_control$tau_shape #Change to 0.0001
  # #d_0 <- prior_control$tau_rate
  
  # beta_lambda<-c(mu_a, rep(1, (K*n_covariates) * (K*n_covariates) + 1) / 2)))
  
  #Regular
  beta_lambda <-c(mu_a, prior_control$mu_log_sig_sq_0, rep(sigma_a, (((K*n_covariates +n_tracers) * ((K*n_covariates +n_tracers) +1))/2)))
  
  #Diag
  #beta_lambda <-c(mu_a, rep(1, K*n_covariates))
  
  # lambda <- c((beta_lambda),
  #   rep(1, n_tracers * 2)
  # )
  lambda <- c(beta_lambda)
  # ,
  #             prior_control$tau_shape, 
  #             prior_control$tau_rate)
  
  
  ll = length(lambda)
  
  lambdares <- c(rep(NA, ll))
  
  
  # thetares <- matrix(rep(NA, ((K * n_covariates + n_tracers * cosimmr_in$n_obs) * n_output)),
  #                    ncol = (K * n_covariates + n_tracers * cosimmr_in$n_obs),
  #                    nrow = n_output
  # )
  
  mylist <- list#vector("list", length = cosimmr_in$n_groups)
  
  # names(mylist) <- levels(cosimmr_in$group)
  
  p_fun <- function(x) exp(x) / sum(exp(x))
  
  # Loop through all the groups
  # for (i in 1:cosimmr_in$n_groups) {
  #  if (cosimmr_in$n_groups > 1) message("\nRunning for group", levels(cosimmr_in$group)[i], "\n\n")
  
  #   curr_rows <- which(cosimmr_in$group_int == i)
  #  curr_mix <- cosimmr_in$mixtures[curr_rows, , drop = FALSE]
  
  # Determine if a single observation or not
  if (nrow(cosimmr_in$mixtures) == 1) {
    message("Only 1 mixture value, performing a simmr solo run...\n")
    solo <- TRUE
    beta_prior = c(rep(100, n_tracers))
  } else {
    solo <- FALSE
    beta_prior = prior_control$tau_rate
  }
  
  n_tracers <- cosimmr_in$n_tracers
  n_sources <- cosimmr_in$n_sources
  s_names <- cosimmr_in$source_names
  K <- cosimmr_in$n_sources
  source_means <- cosimmr_in$source_means
  source_sds <- cosimmr_in$source_sds
  correction_means <- cosimmr_in$correction_means
  correction_sds <- cosimmr_in$correction_sds
  concentration_means <- cosimmr_in$concentration_means
  y <- cosimmr_in$mixtures
  n = nrow(y)
  
  lambdaprior <- c(lambda)
  
  
  #So we take sigma_beta, get the precision, solve, then get transpose so its lower tri, then convert to vector
  # lambdastart <-c(mu_beta_ordered, vec_sig,
  #                 alpha1, beta1, alpha2, beta2)#c(rep(0, length(lambda)))
  #This is only for temporarily using the JAGS starting values
  # lambdastart <- c(1.367, -0.8, -0.373, -0.169, -0.65, -0.189,  0.002, 0.897, -1.482, 0.947, 0.562, 0.195, vec_c1, vec_c2, vec_c3, rep(1, length(lambda - 42)))
  #lambdastart = c(rep(10, length(lambda)))
  # lambdastart = list_1$l_f_j
  lambdastart = lambdaprior
  #lambdastart[91:94] = c(39,50,1,11)
  
  l_l = length(lambdastart)
  #Temporary - want to change the last 4 values to be higher
  # lambdastart[(l_l-3):(l_l)] = lambdastart[(l_l-3):(l_l)] * 10
  lambdastart = lambdaprior * 1
  ### TEMPORARY
  #lambdastart = c(rep(5, length(lambdaprior)))
  
  mu_beta_prior = matrix(prior_control$mu_0, nrow = cosimmr_in$n_covariates)
  sigma_beta_prior = matrix(c(rep(prior_control$sigma_0, cosimmr_in$n_covariates * cosimmr_in$n_sources)), nrow = cosimmr_in$n_covariates)
  
  
  lambdares <- run_VB_cpp(
    lambdastart, K, n_tracers, n_covariates, n, beta_prior, 
    concentration_means,
    source_means, correction_means, correction_sds,
    source_sds, y, (x_scaled), ffvb_control$S,
    ffvb_control$P, ffvb_control$beta_1,
    ffvb_control$beta_2, ffvb_control$tau,
    ffvb_control$eps_0, ffvb_control$t_W,
    prior_control$tau_shape, solo, sigma_beta_prior, mu_beta_prior)
  
  
  iteration = (lambdares$iteration) 
  a = which.max(lambdares$mean_LB_bar[(ffvb_control$t_W +2):iteration]) # +2 here for same reason, start from 0 and then greater than
  use_this = lambdares$lambda_save[(ffvb_control$t_W +1)+a,]
  #use_this = lambdares$lambda_max_LB
  #lambda_from_JAGS = use_this
  #Temporary to check sim_theta
  #use_this = lambdastart
  
  # thetares <- sim_thetacpp0(n_output, lambdares, K, n_tracers, n_covariates, solo)
  kappa2 = rMVNormCpp(K*n_covariates + n_tracers, c(rep(0, n_output)), diag(n_output))
  thetares <- sim_thetacpp(n_output, use_this, K, n_tracers, n_covariates, solo, kappa2)
  
  #p <- t(apply(x_scaled %*% thetares[(1 + n_output * (i - 1)):(n_output * i), 1:K*n_covariates], 1, p_fun))
  
  
  sigma <- sqrt(exp((thetares[,(K*n_covariates + 1):(K*n_covariates + n_tracers)])))
  
  p_sample = array(NA, dim = c(cosimmr_in$n_obs, n_output, K))
  p_mean_sample = matrix(NA, nrow = n_output, ncol = K)
  
  #beta1 = array(thetares[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))
  #The only thing I can potentially think of is that beta is filling wrong
  #So change to a matrix and see if that fixes the problem??
  beta = thetares[,1:(n_covariates * K)]
  
  f <- array(NA, dim = c(cosimmr_in$n_obs, K, n_output)) 
  #f_mean_sample <- array(NA, dim = c(1, K, n_output)) 
  f_mean_sample = matrix(NA, nrow = K, ncol = n_output)
  
  if(cosimmr_in$intercept == TRUE){
    x_sample_mean = c(1, rep(0, (ncol(cosimmr_in$x_scaled) - 1)))
  } else if(cosimmr_in$intercept == FALSE){
    x_sample_mean = c(rep(0, (ncol(cosimmr_in$x_scaled))))
  }
  
  # for(s in 1:n_output){
  #   f[,,s] = as.matrix(x_scaled) %*% beta[s,,]
  # }
  
  #I'm not sure if this works or not
  #We want beta to be the 4 sources for each covariate, 
  #It should be n_covariates * K 
  #n_output here is 3600 cause I keep confusing myself
  for(s in 1:n_output){
    f[,,s] = as.matrix(cosimmr_in$x_scaled) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE)
    # f[,,s] = as.matrix(cosimmr_in$original_x) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE)
    f_mean_sample[,s] = matrix(x_sample_mean, nrow = 1) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE) 
  }
  
  for(j in 1:n_output){
    for (n_obs in 1:cosimmr_in$n_obs) {
      p_sample[n_obs,j, ] <- exp(f[n_obs,1:K, j]) / (sum((exp(f[n_obs,1:K, j]))))
    }
  }
  
  for(j in 1:n_output){
    p_mean_sample[j,] = exp(f_mean_sample[1:K,j]) /(sum(exp(f_mean_sample[1:K,j])))
  }
  
  ########################
  ##vs old way
  # beta1 = array(thetares[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))
  # f1 <- array(NA, dim = c(cosimmr_in$n_obs, K, n_output))
  # for(s in 1:n_output){
  #   f1[,,s] = as.matrix(x_scaled) %*% beta1[s,,]
  # }
  # p_sample1 = array(NA, dim = c(cosimmr_in$n_obs, n_output, K))
  # for(j in 1:n_output){
  #   for (n_obs in 1:cosimmr_in$n_obs) {
  #     p_sample1[n_obs,j, ] <- exp(f1[n_obs,1:K, j]) / (sum((exp(f1[n_obs,1:K, j]))))
  #   }
  # }
  
  
  #Not sure best way to do this??
  #colnames(p[1,,]) <- cosimmr_in$source_names
  
  # sims.matrix <- cbind(
  #   p,
  #   sigma
  # )
  
  # colnames(sims.matrix) <- c(cosimmr_in$source_names, colnames(cosimmr_in$mixtures))
  
  mylist <- list(
    source_names = cosimmr_in$source_names,
    theta = thetares,
    groupnames = cosimmr_in$group,
    lambdares = lambdares,#$lambda,
    # mean_LB = lambdares$mean_LB,
    beta = beta,
    BUGSoutput = list(
      sims.list = list(
        p = p_sample,
        p_mean = p_mean_sample,
        sigma = sigma
      )),
    model = list(data = list(
      mu_f_mean = c(mu_a),
      sigma_f_sd = c(rep(sigma_a, K)),
      ffvb_control = ffvb_control,
      prior_control = prior_control
    ))
  )
  
  
  output_all <- list(input = cosimmr_in, output = mylist)
  
  class(output_all) <- c("cosimmr_output", "ffvb")
  
  return(output_all)}
  else(message("Please select a valid error type"))
}