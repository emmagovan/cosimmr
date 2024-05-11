#'Predicts proportion of each source in a mixture, based on values provided for covariates
#'
#'
#'
#' @param object An object of class \code{cosimmr_output} created via the 
#' function \code{\link{cosimmr_ffvb}}
#' @param x_pred A data.frame of covariate values that the user wishes
#' to predict source proportions for, provided in the same order that the 
#' original covariance matrix was. Important for this to be a data.frame otherwise 
#' numeric values can be set as characters and this causes incorrect calculations.
#' @param n_output the number of posterior samples to generate. Defaults to 3600.
#' @param ... Other arguments (not used)
#'
#' @author Emma Govan <emmagovan@@gmail.com> Andrew Parnell
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
#' @importFrom R2jags jags
#'
#' @examples
#' \donttest{
#' ## See the package vignette for a detailed run through of these 4 examples
#'
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' data(geese_data_day1)
#' cov_1 = c(1,2,3,2,3,1,1,1,2)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ cov_1,
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
#'
#' # Plot
#' plot(simmr_1_out, type = "isospace")
#' plot(simmr_1_out, type = "beta_histogram", cov_name = "cov_1")
#' 
#' x_pred = data.frame(cov_1 = c(1,5))
#'
#'pred_array<-predict(simmr_1_out, x_pred)
#'
#' }
#' @export
predict.cosimmr_output <- function(object,
                          x_pred,
                          n_output = 3600,
                          ...) {
  
#Makes sure the object is the correct class
  if(inherits(object, "cosimmr_output") == TRUE){
    
    ##Need to add something in here to check if its numeric data or categogrical
    
    
    #It has to be a data.frame if theres numeric and categorical data otherwise matrices
    #just turn everything into categorical so it wont work
    if(inherits(x_pred, "data.frame") == FALSE) stop("x_pred must be of type `data.frame` to make predictions with")
    
    scale_x = object$input$scale_x

  
  #Creating a max and min vector so we can check if the new x_pred falls outside the range of the original data
  #Create vectors here but do comparison after all data has been scaled
  max_vec = c(rep(NA, ncol(object$input$x_scaled)))
  min_vec = c(rep(NA, ncol(object$input$x_scaled)))
  
   for(i in 1:(ncol(object$input$x_scaled))){
     max_vec[i] = max(object$input$x_scaled[,i])
     min_vec[i] = min(object$input$x_scaled[,i])
   }
  

  
  
  thetares= object$output$theta
  K = object$input$n_sources
  n_tracers = object$input$n_tracers
  n_covariates = ncol(object$input$x_scaled)
  mixtures = object$input$mixtures
  
  #So now what we want to do is to add x_pred onto the bottom of the original x matrix,
  #scale all, then remove original x mat
  
  original_x = data.frame(object$input$covariates_df)
  
  #Add check that col names match
  #Otherwise next line wont work
  
 # if(colnames(x_pred) != colnames(original_x))stop("Column names of original data and data you wish to predict with do not match. Please fix and rerun.")
  colnames(original_x) = colnames(x_pred)
  
  new_x = rbind(original_x, x_pred)
  
  
  # This adds a column of ones for predicting with if the original had an
  #intercept - want to do this after 
  # if(simmr_out$input$intercept == TRUE){
  #   new_x = cbind(c(rep(1,nrow(new_x))), new_x)
  # } else if(simmr_out$input$intercept == FALSE){
  #   new_x = new_x
  # }
  
  #Now scale
  
  if(nrow(mixtures) == 1){
    #This is if its just 1 entry
    new_x == new_x
  } else{
    if(scale_x == TRUE){
      if(object$input$intercept == TRUE){
        # Original code
        ncol_scaled =  (ncol(stats::model.matrix(~ ., data=new_x))) - 1
        scaled_full_mat = matrix(scale(stats::model.matrix(~ ., data=new_x), 
                           center = c(1,object$input$scaled_center),
                          scale = c(1, object$input$scaled_scale))[,-c(1)], ncol = ncol_scaled)
        scaled_full_mat = cbind(c(rep(1,nrow(scaled_full_mat))), scaled_full_mat)
        
        x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], ncol = ncol(scaled_full_mat))
        
      }else if(object$input$intercept == FALSE){
        scaled_full_mat = scale(stats::model.matrix(~ . -1, data=new_x), 
                                center = object$input$scaled_center,
                                scale = object$input$scaled_scale)
      
        
        x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], ncol = ncol(scaled_full_mat))
        
      }
      
    }else if(scale_x == FALSE){
      if(object$input$intercept == TRUE){
        scaled_full_mat = (stats::model.matrix(~ ., data=new_x))
        
        x_pred_mat = scaled_full_mat[-c(1:nrow(original_x)),]
      }else if(object$input$intercept == FALSE){
        scaled_full_mat = stats::model.matrix(~ .-1, data=new_x)
        
        x_pred_mat = scaled_full_mat[-c(1:nrow(original_x)),]
      }

    }
  }
  
  
  
  #Checks that all the values are above or equal to the min and below or equal to the max
  for(j in 1:(nrow(x_pred_mat))){
    for(i in 1:(ncol(object$input$x_scaled))){
      if(x_pred_mat[j,i] >= min_vec[i] & x_pred_mat[j,i] <= max_vec[i]){
       #message("Data falls within range of data used in original model, okay to predict with")
        print_err = FALSE
      } else(print_err = TRUE)
    }
  }
  
  #This is separate because otherwise its inside the loop and it prints a bunch of times
  if(print_err){message("Please note: The data you wish to predict with falls outside the range of data used in the original model")}
  
  
  
  sigma <- (sqrt(exp(thetares[,(K*n_covariates + 1):(K*n_covariates + n_tracers)])))
  
  p_sample = array(NA, dim =  c(nrow(x_pred_mat), n_output, K))
  
  beta = thetares[,1:(n_covariates * K)]
  
  f <- array(NA, dim = c(nrow(x_pred_mat), K, n_output)) 
  
  for(s in 1:n_output){
    f[,,s] = as.matrix(x_pred_mat) %*% matrix(beta[s,], nrow = n_covariates, ncol = K, byrow = TRUE)
}
  
  for(j in 1:n_output){
    for (n_obs in 1:nrow(x_pred_mat)) {
      p_sample[n_obs,j, ] <- exp(f[n_obs,1:K, j]) / (sum((exp(f[n_obs,1:K, j]))))
    }
  }
  
  out<-list(
    p = p_sample,
    beta = beta,
    sigma = sigma,
    input = object$input,
    theta = thetares
  )
  
  class(out) = "cosimmr_pred_out"
 return(out)
  } else (message("Can only predict using cosimmr_output object generated from
cosimmr_ffvb function"))
}
