
#' Run a \code{cosimmr_input} object through the Fixed Form Variational
#' Bayes(FFVB) function
#'
#' This is the main function of cosimmr. It takes a \code{cosimmr_input} object
#' created via \code{\link{cosimmr_load}}, runs it in fixed form
#' Variational Bayes to determine the dietary proportions, and then
#' outputs a \code{cosimmr_output} object for further analysis and plotting
#' via \code{\link{plot.cosimmr_output}}.
#'

#' @param simmr_in An object created via the function \code{\link{cosimmr_load}}
#' @param prior_control A list of values including arguments named \code{mu_0}
#' (prior for mu), and \code{sigma_0} (prior for sigma).
#' @param ffvb_control A list of values including arguments named \code{n_output}
#' (number of rows in theta output), \code{S} (number of samples taken at each
#' iteration of the algorithm), \code{P} (patience parameter), \code{beta_1}
#' and \code{beta_2} (adaptive learning weights), \code{tau} (threshold for
#' exploring learning space), \code{eps_0} (fixed learning rate),
#' \code{t_W} (rolling window size)
#' @return An object of class \code{cosimmr_output} with two named top-level
#' components: \item{input }{The \code{cosimmr_input} object given to the
#' \code{simmr_ffvb} function} \item{output }{A set of outputs produced by
#' the FFVB function. These can be analysed using the
#' \code{\link{summary.cosimmr_output}} and \code{\link{plot.cosimmr_output}}
#' functions.}
#'
#' @author Emma Govan <emma.govan.2021@@mumail.ie>
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
#' \dontrun{
#' ## See the package vignette for a detailed run through of these 4 examples
#'
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' data(geese_data_day1)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ c(1,2,3,2,1,3,2,1,2),
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
#' plot(simmr_1)
#'
#' # Print
#' simmr_1
#'
#' # FFVB run
#' simmr_1_out <- cosimmr_ffvb(simmr_1)
#'
#' # Print it
#' print(simmr_1_out)
#'
#' # Summary
#' summary(simmr_1_out, type = "correlations")
#' summary(simmr_1_out, type = "statistics")
#' ans <- summary(simmr_1_out, type = c("quantiles", "statistics"))
#'
#' # Plot
#' plot(simmr_1_out, type = "beta_boxplot")
#' plot(simmr_1_out, type = "beta_histogram")
#'
#' # Compare two sources
#' compare_sources(simmr_1_out, source_names = c("Zostera", "Enteromorpha"))
#'
#' # Compare multiple sources
#' compare_sources(simmr_1_out)
#'
#' #####################################################################################
#'
#' # A version with just one observation
#' data(geese_data_day1)
#' simmr_2 <- with(
#'   geese_data_day1,
#'   simmr_loadcov(
#'     formula = mixtures[1, , drop = FALSE] ~ c(1),
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
#' plot(simmr_2)
#'
#' # FFVB run - automatically detects the single observation
#' simmr_2_out <- simmr_ffvb(simmr_2)
#'
#' # Print it
#' print(simmr_2_out)
#'
#' # Summary
#' summary(simmr_2_out)
#' ans <- summary(simmr_2_out, type = c("quantiles"))
#'
#' # Plot
#' plot(simmr_2_out)
#' plot(simmr_2_out, type = "boxplot")
#' plot(simmr_2_out, type = "histogram")
#' plot(simmr_2_out, type = "density")
#' plot(simmr_2_out, type = "matrix")
#'
#' #####################################################################################
#'
#' # Data set 2: 3 isotopes (d13C, d15N and d34S), 30 observations, 4 sources
#' data(simmr_data_2)
#' simmr_3 <- with(
#'   simmr_data_2,
#'   cosimmr_load(
#'     formula = mixtures ~ c(rep(1, ncol(mixtures))),
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means
#'   )
#' )
#'
#' # Get summary
#' print(simmr_3)
#'
#' # Plot 3 times
#' plot(simmr_3)
#' plot(simmr_3, tracers = c(2, 3))
#' plot(simmr_3, tracers = c(1, 3))
#' # See vignette('simmr') for fancier axis labels
#'
#' # FFVB run
#' simmr_3_out <- cosimmr_ffvb(simmr_3)
#'
#' # Print it
#' print(simmr_3_out)
#'
#' # Summary
#' summary(simmr_3_out)
#' summary(simmr_3_out, type = "quantiles")
#' summary(simmr_3_out, type = "correlations")
#'
#' # Plot
#' plot(simmr_3_out)
#' plot(simmr_3_out, type = "boxplot")
#' plot(simmr_3_out, type = "histogram")
#' plot(simmr_3_out, type = "density")
#' plot(simmr_3_out, type = "matrix")
#'
#' ################################################################
#'}
#' @export cosimmr_ffvb
cosimmr_ffvb <- function(simmr_in,
                         prior_control = list(
                           mu_0 = rep(0, simmr_in$n_sources * simmr_in$n_covariates),
                           sigma_0 = 1,
                           tau_shape = rep(1, simmr_in$n_tracers),
                           tau_rate = rep(1, simmr_in$n_tracers)
                         ),
                         ffvb_control = list(
                           n_output = 3600,
                           S = 1000,
                           P = 50,
                           beta_1 = 0.9,
                           beta_2 = 0.9,
                           tau = 100,
                           eps_0 = 0.005,
                           t_W = 100
                         )) {
  # Throw a warning if less than 4 observations in a group - 1 is ok as it wil do a solo run
  #  if (min(table(simmr_in$group)) > 1 & min(table(simmr_in$group)) < 4) warning("At least 1 group has less than 4 observations - either put each observation in an individual group or use informative prior information")
  
  output <- list #vector("list", length = simmr_in$n_groups)
  # names(output) <- levels(simmr_in$group)
  K <- simmr_in$n_sources
  n_tracers <- simmr_in$n_tracers
  n_output <- ffvb_control$n_output
  mu_a <- prior_control$mu_0
  sigma_a <- prior_control$sigma_0
  n_covariates<- simmr_in$n_covariates
  x_scaled = simmr_in$x_scaled
  
  c_0 <- prior_control$tau_shape #Change to 0.0001
  # #d_0 <- prior_control$tau_rate
  
 # beta_lambda<-c(mu_a, rep(1, (K*n_covariates) * (K*n_covariates) + 1) / 2)))

  #Regular
beta_lambda <-c(mu_a, rep(1,(((K*n_covariates) * (K*n_covariates +1))/2)))

  #Diag
#beta_lambda <-c(mu_a, rep(1, K*n_covariates))
  
  lambda <- c((beta_lambda),
    rep(1, n_tracers * 2)
  )
  
  ll = length(lambda)
  
  lambdares <- c(rep(NA, ll))
  
  thetares <- matrix(rep(NA, ((K * n_covariates + n_tracers) * n_output)),
                     ncol = (K * n_covariates + n_tracers),
                     nrow = n_output
  )
  
  mylist <- list#vector("list", length = simmr_in$n_groups)
  
  # names(mylist) <- levels(simmr_in$group)
  
  p_fun <- function(x) exp(x) / sum(exp(x))
  
  # Loop through all the groups
  # for (i in 1:simmr_in$n_groups) {
  #  if (simmr_in$n_groups > 1) message("\nRunning for group", levels(simmr_in$group)[i], "\n\n")
  
  #   curr_rows <- which(simmr_in$group_int == i)
  #  curr_mix <- simmr_in$mixtures[curr_rows, , drop = FALSE]
  
  # Determine if a single observation or not
  if (nrow(simmr_in$mixtures) == 1) {
    message("Only 1 mixture value, performing a simmr solo run...\n")
    solo <- TRUE
    beta_prior = c(rep(100, n_tracers))
  } else {
    solo <- FALSE
    beta_prior = prior_control$tau_rate
  }
  
  n_tracers <- simmr_in$n_tracers
  n_sources <- simmr_in$n_sources
  s_names <- simmr_in$source_names
  K <- simmr_in$n_sources
  source_means <- simmr_in$source_means
  source_sds <- simmr_in$source_sds
  correction_means <- simmr_in$correction_means
  correction_sds <- simmr_in$correction_sds
  concentration_means <- simmr_in$concentration_means
  y <- simmr_in$mixtures
  n = nrow(y)
  
  lambdastart <- c(lambda)
  
  lambdares <- run_VB_cpp(
    lambdastart, K, n_tracers, n_covariates, n, beta_prior, concentration_means,
    source_means, correction_means, correction_sds,
    source_sds, y, (x_scaled), ffvb_control$S,
    ffvb_control$P, ffvb_control$beta_1,
    ffvb_control$beta_2, ffvb_control$tau,
    ffvb_control$eps_0, ffvb_control$t_W,
    c_0, solo
  )
  
  
  
  #thetares <- sim_thetacpp(n_output, lambdares, K, n_tracers, n_covariates, solo)
  thetares <- sim_thetacpp(n_output, lambdares$lambda, K, n_tracers, n_covariates, solo)
  
  #p <- t(apply(x_scaled %*% thetares[(1 + n_output * (i - 1)):(n_output * i), 1:K*n_covariates], 1, p_fun))
  
  
  sigma <- (1/sqrt(thetares[,(K*n_covariates + 1):(K*n_covariates + n_tracers)]))
  
  p_sample = array(NA, dim = c(simmr_in$n_obs, n_output, K))
  p_mean_sample = matrix(NA, nrow = n_output, ncol = K)
  
  #beta1 = array(thetares[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))
  #The only thing I can potentially think of is that beta is filling wrong
  #So change to a matrix and see if that fixes the problem??
  beta = thetares[,1:(n_covariates * K)]
  
  f <- array(NA, dim = c(simmr_in$n_obs, K, n_output)) 
  #f_mean_sample <- array(NA, dim = c(1, K, n_output)) 
  f_mean_sample = matrix(NA, nrow = K, ncol = n_output)
  
  if(simmr_in$intercept == TRUE){
    x_sample_mean = c(1, rep(0, (ncol(simmr_in$x_scaled) - 1)))
  } else if(simmr_in$intercept == FALSE){
    x_sample_mean = c(rep(0, (ncol(simmr_in$x_scaled))))
  }
  
  # for(s in 1:n_output){
  #   f[,,s] = as.matrix(x_scaled) %*% beta[s,,]
  # }
  
  #I'm not sure if this works or not
  #We want beta to be the 4 sources for each covariate, 
  #It should be n_covariates * K 
  #n_output here is 3600 cause I keep confusing myself
  for(s in 1:n_output){
    
    f[,,s] = as.matrix(x_scaled) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE)
    f_mean_sample[,s] = matrix(x_sample_mean, nrow = 1) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE) 
  }
  
  for(j in 1:n_output){
    for (n_obs in 1:simmr_in$n_obs) {
      p_sample[n_obs,j, ] <- exp(f[n_obs,1:K, j]) / (sum((exp(f[n_obs,1:K, j]))))
    }
  }
  
  for(j in 1:n_output){
    p_mean_sample[j,] = exp(f_mean_sample[1:K,j]) /(sum(exp(f_mean_sample[1:K,j])))
  }
  
  ########################
  ##vs old way
  # beta1 = array(thetares[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))
  # f1 <- array(NA, dim = c(simmr_in$n_obs, K, n_output))
  # for(s in 1:n_output){
  #   f1[,,s] = as.matrix(x_scaled) %*% beta1[s,,]
  # }
  # p_sample1 = array(NA, dim = c(simmr_in$n_obs, n_output, K))
  # for(j in 1:n_output){
  #   for (n_obs in 1:simmr_in$n_obs) {
  #     p_sample1[n_obs,j, ] <- exp(f1[n_obs,1:K, j]) / (sum((exp(f1[n_obs,1:K, j]))))
  #   }
  # }
  
  
  #Not sure best way to do this??
  #colnames(p[1,,]) <- simmr_in$source_names
  
  # sims.matrix <- cbind(
  #   p,
  #   sigma
  # )
  
  # colnames(sims.matrix) <- c(simmr_in$source_names, colnames(simmr_in$mixtures))
  
  mylist <- list(
    source_names = simmr_in$source_names,
    theta = thetares,
    groupnames = simmr_in$group,
    lambdares = lambdares$lambda,
    mean_LB = lambdares$mean_LB,
    beta = beta,
    BUGSoutput = list(
      sims.list = list(
        p = p_sample,
        p_mean = p_mean_sample,
        sigma = sigma
      )),
    model = list(data = list(
      mu_f_mean = c(mu_a),
      sigma_f_sd = c(rep(sigma_a, K))
    ))
  )
  
  
  output_all <- list(input = simmr_in, output = mylist)
  
  class(output_all) <- c("cosimmr_output", "ffvb")
  
  return(output_all)
}
