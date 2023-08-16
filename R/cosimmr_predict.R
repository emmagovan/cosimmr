#'Predicts proportion of each source in a mixture, based on values provided for covariates
#'
#'
#'
#' @param simmr_out An object created via the function \code{\link{cosimmr_ffvb}}
#' @param x_pred A vector or matrix of covariate values that the user wishes
#' to predict source proportions for, provided in the same order that the 
#' original covariance matrix was

#'
#' @author Emma Govan <emma.govan.2021@mumail.ie>
#'
#' @seealso \code{\link{cosimmr_load}} for creating objects suitable for this
#' function,  and
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
#'     mixtures = mixtures,
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
#' plot(simmr_1_out, type = "boxplot")
#' plot(simmr_1_out, type = "histogram")
#' plot(simmr_1_out, type = "density")
#' plot(simmr_1_out, type = "matrix")
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
#'   cosimmr_load(
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
#' simmr_2_out <- cosimmr_ffvb(simmr_2)
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
#' }
#' @export
simmrcov_predict <- function(simmr_out,
                          x_pred,
                          n_output = 3600) {
  

  
 # x_pred_mat = matrix(x_pred, ncol = n_covariates)
  
  #Check x_pred has same number of columns as original x
  #This means they have to input x_pred as a matrix
  #Otherwise this wont work
  x_pred_mat = x_pred 
  
if(ncol(x_pred_mat) != ncol(simmr_out$input$x_scaled)) stop("The matrix of values you wish to make predictions for does not have the same number of entries as the original covariance matrix. Please fix and rerun.")

  #Not sure if we want to include this - do we want to ensure the columns are named??
#if(colnames(x_pred_mat) != colnames(simmr_out$input$x_scaled)) stop("The column names for the original covariates and the new values you wish to make predictions for are not the same")
  
  
  #Creating a max and min vector so we can check if the new x_pred falls outside the range of the original data
  #Create vectors here but do comparison after all data has been scaled
  max_vec = c(rep(NA, ncol(simmr_out$input$x_scaled)))
  min_vec = c(rep(NA, ncol(simmr_out$input$x_scaled)))
  
  for(i in 1:(ncol(simmr_out$input$x_scaled))){
    max_vec[i] = max(simmr_out$input$x_scaled[,i])
    min_vec[i] = min(simmr_out$input$x_scaled[,i])
  }
  

  
  
  thetares= simmr_out$output$theta
  K = simmr_out$input$n_sources
  n_tracers = simmr_out$input$n_tracers
  n_covariates = ncol(simmr_out$input$x_scaled)
  
  #This makes sure that the new inputs are scaled in the same way the old ones
  #were - so if they are scaled there they're scaled here using the same
  #scaling values (im sure theres a technical term for this but I don't know it)
  if(simmr_out$input$scale_x == FALSE){
  x_pred_mat = matrix(x_pred, ncol = n_covariates)
  } else if(simmr_out$input$scale_x == TRUE){
    if(simmr_out$input$intercept == TRUE){
      x_pred_mat = cbind(x_pred_mat[,1], scale(x_pred_mat[,2:ncol(x_pred_mat)], 
                                               center = simmr_out$input$scaled_center,
                                               scale = simmr_out$input$scaled_scale))
    }
    else if(simmr_out$input$intercept == FALSE){
      x_pred_mat = scale(x_pred_mat, 
                         center = simmr_out$input$scaled_center,
                         scale = simmr_out$input$scaled_scale)
    }
  }
  
  #Checks that all the values are above or equal to the min and below or equal to the max
  for(j in 1:(nrow(x_pred_mat))){
    for(i in 1:(ncol(simmr_out$input$x_scaled))){
      if(x_pred_mat[j,i] >= min_vec[i] & x_pred_mat[j,i] <= max_vec[i]){
       #message("Data falls within range of data used in original model, okay to predict with")
        print_err = FALSE
      } else(print_err = TRUE)
    }
  }
  
  #This is separate because otherwise its inside the loop and it prints a bunch of times
  if(print_err == TRUE){message("Please note: The data you wish to predict with falls outside the range of data used in the original model")}
  
  
  
  sigma <- (1/sqrt(thetares[,(K*n_covariates + 1):(K*n_covariates + n_tracers)]))
  
  p_sample = array(NA, dim =  c(nrow(x_pred_mat), n_output, K))
  
  beta = array(thetares[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))
  
  f <- array(NA, dim = c(nrow(x_pred_mat), K, n_output)) 
  
  for(s in 1:n_output){
    f[,,s] = as.matrix(x_pred_mat) %*% beta[s,,]
  }
  
  for(j in 1:n_output){
    for (n_obs in 1:nrow(x_pred_mat)) {
      p_sample[n_obs,j, ] <- exp(f[n_obs,1:K, j]) / (sum((exp(f[n_obs,1:K, j]))))
    }
  }
 return(p_sample)
}
