#' Function to load in cosimmr data and check for errors
#'
#' This function takes in the mixture data, food source means and standard
#' deviations, and (optionally) correction factor means and standard
#' deviations, and concentration proportions. It performs some (non-exhaustive)
#' checking of the data to make sure it will run through simmr. It outputs an
#' object of class \code{cosimmr_input}.
#'
#' For standard stable isotope mixture modelling, the mixture matrix will
#' contain a row for each individual and a column for each isotopic value.
#' \code{cosimmr} will allow for any number of isotopes and any number of
#' observations, within computational limits. The source means/sds should be
#' provided for each food source on each isotope. The correction means (usually
#' trophic enrichment factors) can be set as zero if required, and should be of
#' the same shape as the source values. The concentration dependence means
#' should be estimated values of the proportion of each element in the food
#' source in question and should be given in proportion format between 0 and 1.
#' At present there is no means to include concentration standard deviations.
#'
#' @param formula Formula giving in form y ~ x where y is a vector or matrix
#' of mixture values and x is a vector or matrix of covariates
#' @param source_names The names of the sources given as a character string
#' @param source_means The means of the source values, given as a matrix where
#' the number of rows is the number of sources and the number of columns is the
#' number of tracers
#' @param source_sds The standard deviations of the source values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers
#' @param correction_means The means of the correction values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. If not provided these are set to 0.
#' @param correction_sds The standard deviations of the correction values,
#' given as a matrix where the number of rows is the number of sources and the
#' number of columns is the number of tracers. If not provided these are set to
#' 0.
#' @param concentration_means The means of the concentration values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. These should be between 0 and 1. If not
#' provided these are all set to 1.
#' @param scale_x Whether or not you wish to scale the x values provided, or run 
#' the model using the original x values. Defaults to TRUE.
#'
#' @import checkmate
#'
#'
#' @return An object of class \code{cosimmr_input} with the following elements:
#' \item{mixtures }{The mixture data} \item{source_names }{Source means}
#' \item{sources_sds }{Source standard deviations} \item{correction_means
#' }{Correction means} \item{correction_sds }{Correction standard deviations}
#' \item{concentration_means }{Concentration dependence means} \item{n_obs
#' }{The number of observations} \item{n_tracers }{The number of
#' tracers/isotopes} \item{n_sources }{The number of sources} \item{n_groups
#' }{The number of groups}
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
#' @seealso See \code{\link{cosimmr_ffvb}} for complete examples.
#' @examples
#' \dontrun{
#'
#' # A simple example with 10 observations, 2 tracers and 4 sources
#' data(geese_data_day1)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ 1,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means,
#'     scale_x = TRUE
#'   )
#' )
#' 
#' 
#'
#' print(simmr_1)
#' }
#' @export cosimmr_load
cosimmr_load <- function(formula,
source_names,
source_means,
source_sds,
correction_means = NULL,
correction_sds = NULL,
concentration_means = NULL,
scale_x = TRUE) {
  # Function to load in data for simmr and check whether it's appropriate for running through simmr_mcmc
  
  
  
  #Possibly need covariate data frame
  #c_df = data.frame(colour_cat = colour_cat, letter_cat = letter_cat, numeric_cov = numeric_cov)
  # Go through each object and check that it matches the requirements
  mixtures = as.matrix(stats::model.frame(formula)[,1])
  original_x = stats::model.matrix(formula)
  
  #Need to add some method here to keep the column names
  covariates = (stats::model.frame(formula)[,-1])
  
  
  # model.matrix(~ ., data=c_df, 
  # contrasts.arg = lapply(c_df[sapply(c_df, is.factor)],
  #                        contrasts,
  #                        contrasts=FALSE))
  # 
  
  if(nrow(mixtures) == 1){
    #This is if its just 1 entry
    x_scaled = stats::model.matrix(formula)
  } 
  if(scale_x == TRUE){
    if(stats::sd(stats::model.matrix(formula)[,1]) == 0){
      if(ncol(stats::model.matrix(formula)) == 1){
        intercept = TRUE
        scaled_mat = (stats::model.matrix(formula))
        x_scaled = scaled_mat
        colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))
        
        scaled_center = attr(scaled_mat, "scaled:center")
        
        scaled_scale = attr(scaled_mat, "scaled:scale")
        print("Cannot scale when only using row of 1s")
      } else if(ncol(stats::model.matrix(formula)) != 1){
        # Original code
        intercept = TRUE
        scaled_mat = scale(stats::model.matrix(formula)[,(2:ncol(stats::model.matrix(formula)))])
        x_scaled = cbind(stats::model.matrix(formula)[,1],
                         scaled_mat)
        colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))
        
        scaled_center = attr(scaled_mat, "scaled:center")
        
        scaled_scale = attr(scaled_mat, "scaled:scale")
        
      }
      
    }else if(stats::sd(stats::model.matrix(formula)[,1]) != 0){
      intercept = FALSE
      scaled_mat = scale(stats::model.matrix(formula))
      x_scaled = scaled_mat
      colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))
      
      scaled_center = attr(scaled_mat, "scaled:center")
      
      scaled_scale = attr(scaled_mat, "scaled:scale")
      
    }
  } else if(scale_x == FALSE){
    
    if(stats::sd(stats::model.matrix(formula)[,1]) == 0){
      intercept = TRUE}else(intercept = FALSE)
    x_scaled = stats::model.matrix(formula)
    scaled_center = NULL
    
    scaled_scale = NULL
    
    
  }
  
  # Write a function that generically tests for any 2D numeric data shape such as matrix, data frame or tibble
  assert_2D_numeric <- function(x,
                                nrows = NULL,
                                ncols = NULL,
                                null.ok = FALSE) {
    assert(
      test_data_frame(x,
                      types = c("double", "numeric"),
                      nrows = nrows,
                      ncols = ncols,
                      null.ok = null.ok
      ),
      test_matrix(x,
                  mode = "numeric",
                  nrows = nrows,
                  ncols = ncols,
                  null.ok = null.ok
      ),
      test_tibble(x,
                  types = c("double", "numeric"),
                  nrows = nrows,
                  ncols = ncols,
                  null.ok = null.ok
      )
    )
  }
  
  # Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
  # assert_matrix(mixtures)
  assert_2D_numeric(mixtures)
  n_obs <- nrow(mixtures)
  n_tracers <- ncol(mixtures)
  
  # Add column names if they're not there
  if (is.null(colnames(mixtures))) {
    colnames(mixtures) <- paste0("tracer", 1:n_tracers)
  }
  
  if (is.null(colnames(x_scaled))) {
    colnames(x_scaled) <- paste0("covariate", 1:n_tracers)
  }
  
  # source_names must be a character vector - the length of it is the number of sources
  assert_character(source_names)
  n_sources <- length(source_names)
  
  # source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
  assert_2D_numeric(source_means,
                    nrows = n_sources,
                    ncols = n_tracers
  )
  # assert_matrix(source_means, nrows = n_sources, ncols = n_tracers)
  assert_2D_numeric(source_sds,
                    nrows = n_sources,
                    ncols = n_tracers
  )
  # assert_matrix(source_sds, nrows = n_sources, ncols = n_tracers)
  assert_2D_numeric(correction_means,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = ifelse(is.null(correction_sds),
                                     TRUE, FALSE
                    )
  )
  # assert_matrix(correction_means,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_sds),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(correction_sds,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = ifelse(is.null(correction_sds),
                                     TRUE, FALSE
                    )
  )
  # assert_matrix(correction_sds,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_means),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(concentration_means,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = TRUE
  )
  # assert_matrix(concentration_means,
  #   nrows = n_sources,
  #   ncols = n_tracers, null.ok = TRUE
  # )
  
  # Fill in correction means
  if (is.null(correction_means)) {
    correction_means <- matrix(0, ncol = n_tracers, nrow = n_sources)
    correction_sds <- matrix(0, ncol = n_tracers, nrow = n_sources)
  }
  
  # concentration_means must be a matrix where all elements are less than 1
  if (is.null(concentration_means)) {
    concentration_means <- matrix(1, ncol = n_tracers, nrow = n_sources)
  } else {
    assert_true(all(concentration_means < 1) & all(concentration_means > 0))
  }
  
  # Check the groups are the right length and structure if given
  
  
  
  # Prepare output and give class
  out <- list(
    mixtures = mixtures,
    x_scaled = x_scaled,
    source_names = source_names,
    source_means = source_means,
    source_sds = source_sds,
    correction_means = correction_means,
    correction_sds = correction_sds,
    concentration_means = concentration_means,
    n_obs = n_obs,
    n_tracers = n_tracers,
    n_sources = n_sources,
    scale_x = scale_x,
    scaled_center = scaled_center,
    scaled_scale = scaled_scale,
    intercept = intercept,
    covariates_df = covariates,
    n_covariates =  ncol(x_scaled),
    original_x = original_x
  )
  
  # Look through to see whether there are any missing values in anything
  if (any(unlist(lapply(out, "is.na")))) {
    warning("Missing values provided for some values. Check your inputs")
  }
  
  class(out) <- "cosimmr_input"
  
  return(out)
}

# cosimmr_load <- function(formula,
#                          source_names,
#                          source_means,
#                          source_sds,
#                          correction_means = NULL,
#                          correction_sds = NULL,
#                          concentration_means = NULL,
#                          scale_x = TRUE) {
#   # Function to load in data for simmr and check whether it's appropriate for running through simmr_mcmc
#   
#   
#   
#   #Possibly need covariate data frame
#   #c_df = data.frame(colour_cat = colour_cat, letter_cat = letter_cat, numeric_cov = numeric_cov)
#   # Go through each object and check that it matches the requirements
#   mixtures = as.matrix(stats::model.frame(formula)[,1])
#   original_x = stats::model.matrix(formula)
#   
#   #Need to add some method here to keep the column names
#   covariates = (stats::model.frame(formula)[,-1])
#   
#   
#   # model.matrix(~ ., data=c_df, 
#   # contrasts.arg = lapply(c_df[sapply(c_df, is.factor)],
#   #                        contrasts,
#   #                        contrasts=FALSE))
#   # 
#   
#   if(nrow(mixtures) == 1){
#     #This is if its just 1 entry
#     x_scaled = stats::model.matrix(formula)
#   } 
#   if(scale_x == TRUE){
#     if(stats::sd(stats::model.matrix(formula)[,1]) == 0){
#       if(ncol(stats::model.matrix(formula)) == 1){
#         intercept = TRUE
#         scaled_mat = (stats::model.matrix(formula))
#         x_scaled = scaled_mat
#         colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))
#         
#         scaled_center = attr(scaled_mat, "scaled:center")
#         
#         scaled_scale = attr(scaled_mat, "scaled:scale")
#         print("Cannot scale when only using row of 1s")
#       } else if(ncol(stats::model.matrix(formula)) != 1){
#         # Original code
#         intercept = TRUE
#         scaled_mat = scale(stats::model.matrix(formula)[,(2:ncol(stats::model.matrix(formula)))])
#         x_scaled = cbind(stats::model.matrix(formula)[,1],
#                          scaled_mat)
#         colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))
#         
#         scaled_center = attr(scaled_mat, "scaled:center")
#         
#         scaled_scale = attr(scaled_mat, "scaled:scale")
#         
#       }
#       
#     }else if(stats::sd(stats::model.matrix(formula)[,1]) != 0){
#       intercept = FALSE
#       scaled_mat = scale(stats::model.matrix(formula))
#       x_scaled = scaled_mat
#       colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))
#       
#       scaled_center = attr(scaled_mat, "scaled:center")
#       
#       scaled_scale = attr(scaled_mat, "scaled:scale")
#       
#     }
#   } else if(scale_x == FALSE){
#     
#     if(stats::sd(stats::model.matrix(formula)[,1]) == 0){
#       intercept = TRUE}else(intercept = FALSE)
#     x_scaled = stats::model.matrix(formula)
#     scaled_center = NULL
#     
#     scaled_scale = NULL
#     
#     
#   }
#   
#   # Write a function that generically tests for any 2D numeric data shape such as matrix, data frame or tibble
#   assert_2D_numeric <- function(x,
#                                 nrows = NULL,
#                                 ncols = NULL,
#                                 null.ok = FALSE) {
#     assert(
#       test_data_frame(x,
#                       types = c("double", "numeric"),
#                       nrows = nrows,
#                       ncols = ncols,
#                       null.ok = null.ok
#       ),
#       test_matrix(x,
#                   mode = "numeric",
#                   nrows = nrows,
#                   ncols = ncols,
#                   null.ok = null.ok
#       ),
#       test_tibble(x,
#                   types = c("double", "numeric"),
#                   nrows = nrows,
#                   ncols = ncols,
#                   null.ok = null.ok
#       )
#     )
#   }
#   
#   # Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
#   # assert_matrix(mixtures)
#   assert_2D_numeric(mixtures)
#   n_obs <- nrow(mixtures)
#   n_tracers <- ncol(mixtures)
#   
#   # Add column names if they're not there
#   if (is.null(colnames(mixtures))) {
#     colnames(mixtures) <- paste0("tracer", 1:n_tracers)
#   }
#   
#   if (is.null(colnames(x_scaled))) {
#     colnames(x_scaled) <- paste0("covariate", 1:n_tracers)
#   }
#   
#   # source_names must be a character vector - the length of it is the number of sources
#   assert_character(source_names)
#   n_sources <- length(source_names)
#   
#   # source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
#   assert_2D_numeric(source_means,
#                     nrows = n_sources,
#                     ncols = n_tracers
#   )
#   # assert_matrix(source_means, nrows = n_sources, ncols = n_tracers)
#   assert_2D_numeric(source_sds,
#                     nrows = n_sources,
#                     ncols = n_tracers
#   )
#   # assert_matrix(source_sds, nrows = n_sources, ncols = n_tracers)
#   assert_2D_numeric(correction_means,
#                     nrows = n_sources,
#                     ncols = n_tracers,
#                     null.ok = ifelse(is.null(correction_sds),
#                                      TRUE, FALSE
#                     )
#   )
#   # assert_matrix(correction_means,
#   #   nrows = n_sources,
#   #   ncols = n_tracers,
#   #   null.ok = ifelse(is.null(correction_sds),
#   #     TRUE, FALSE
#   #   )
#   # )
#   assert_2D_numeric(correction_sds,
#                     nrows = n_sources,
#                     ncols = n_tracers,
#                     null.ok = ifelse(is.null(correction_sds),
#                                      TRUE, FALSE
#                     )
#   )
#   # assert_matrix(correction_sds,
#   #   nrows = n_sources,
#   #   ncols = n_tracers,
#   #   null.ok = ifelse(is.null(correction_means),
#   #     TRUE, FALSE
#   #   )
#   # )
#   assert_2D_numeric(concentration_means,
#                     nrows = n_sources,
#                     ncols = n_tracers,
#                     null.ok = TRUE
#   )
#   # assert_matrix(concentration_means,
#   #   nrows = n_sources,
#   #   ncols = n_tracers, null.ok = TRUE
#   # )
#   
#   # Fill in correction means
#   if (is.null(correction_means)) {
#     correction_means <- matrix(0, ncol = n_tracers, nrow = n_sources)
#     correction_sds <- matrix(0, ncol = n_tracers, nrow = n_sources)
#   }
#   
#   # concentration_means must be a matrix where all elements are less than 1
#   if (is.null(concentration_means)) {
#     concentration_means <- matrix(1, ncol = n_tracers, nrow = n_sources)
#   } else {
#     assert_true(all(concentration_means < 1) & all(concentration_means > 0))
#   }
#   
#   # Check the groups are the right length and structure if given
#   
#   
#   
#   # Prepare output and give class
#   out <- list(
#     mixtures = mixtures,
#     x_scaled = x_scaled,
#     source_names = source_names,
#     source_means = source_means,
#     source_sds = source_sds,
#     correction_means = correction_means,
#     correction_sds = correction_sds,
#     concentration_means = concentration_means,
#     n_obs = n_obs,
#     n_tracers = n_tracers,
#     n_sources = n_sources,
#     scale_x = scale_x,
#     scaled_center = scaled_center,
#     scaled_scale = scaled_scale,
#     intercept = intercept,
#     covariates_df = covariates,
#     n_covariates =  ncol(x_scaled),
#     original_x = original_x
#   )
#   
#   # Look through to see whether there are any missing values in anything
#   if (any(unlist(lapply(out, "is.na")))) {
#     warning("Missing values provided for some values. Check your inputs")
#   }
#   
#   class(out) <- "cosimmr_input"
#   
#   return(out)
# }