#' Summarises the output created with \code{\link{cosimmr_ffvb}}
#'
#' Produces textual summaries and convergence diagnostics for an object created
#' with  \code{\link{cosimmr_ffvb}}. The different
#' options are: 'quantiles' which produces credible intervals
#' for the parameters, 'statistics' which produces means and standard
#' deviations, and 'correlations' which produces correlations between the
#' parameters.
#'
#' The quantile output allows easy calculation of 95 per cent credible
#' intervals of the posterior dietary proportions. The correlations allow the 
#' user to judge which sources are non-identifiable.
#'
#'
#' @param object An object of class \code{cosimmr_output} produced by the
#' function \code{\link{cosimmr_ffvb}}
#' @param type The type of output required. At least none of quantiles', 
#' 'statistics', or 'correlations'.
#' @param obs The observation to generate a summary for. Defaults to 1.
#' @param ...  Not used
#' @return A list containing the following components: 
#' \item{quantiles }{The quantiles of each parameter from the posterior 
#' distribution} \item{statistics }{The means and standard
#' deviations of each parameter} \item{correlations }{The posterior
#' correlations between the parameters} Note that this object is reported
#' silently so will be discarded unless the function is called with an object
#' as in the example below.
#' @author Emma Govan <emmagovan@@gmail.com> Andrew Parnell
#' @seealso See \code{\link{cosimmr_ffvb}}for creating objects suitable for 
#' this function, and many more examples.
#' See also \code{\link{cosimmr_load}} for creating cosimmr objects,
#' \code{\link{plot.cosimmr_input}} for creating isospace plots,
#' \code{\link{plot.cosimmr_output}} for plotting output.
#'
#' @importFrom stats sd cor
#'
#' @examples
#' \donttest{
#' # A simple example with 10 observations, 2 tracers and 4 sources
#'
#' # The data
#' data(geese_data_day1)
#' cosimmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ 1,
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
#'
#' # FFVB run
#' cosimmr_1_out <- cosimmr_ffvb(cosimmr_1)
#'
#' # Summarise
#' summary(cosimmr_1_out) # This outputs all the summaries
#' summary(cosimmr_1_out, type = "quantiles") # Just the diagnostics
#' # Store the output in an
#' ans <- summary(cosimmr_1_out,
#'   type = c("quantiles", "statistics")
#' )
#' }
#' @export
summary.cosimmr_output <-
function(object, type = c("quantiles", "statistics", "correlations"), obs = 1, ...) {
  if (inherits(object, "cosimmr_output") == TRUE) {
    if (inherits(object, "ffvb") == TRUE) {
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)
      
      
      
      # Set up containers
      out_bgr <- out_quantiles <- out_statistics <- out_cor <- vector("list", length = length(obs))
      names(out_bgr) <- paste0("obs_", obs)
      names(out_quantiles) <- paste0("obs_", obs)
      names(out_statistics) <- paste0("obs_", obs)
      names(out_cor) <- paste0("obs_", obs)
      
      # Loop through groups
      for (i in 1:length(obs)) {
        message("\nSummary for Observation ",  obs[i], "\n")
        out_all <- cbind(object$output$BUGSoutput$sims.list$p[obs[i],,],
                         object$output$BUGSoutput$sims.list$sigma)
        colnames(out_all) = c(object$input$source_names, paste0("sd_", colnames(object$input$mixtures)))
        
        # Get cosimmr_outs
        out_quantiles[[i]] <- t(apply(out_all, 2, "quantile", probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
        #  coda:::summary.mcmc.list(cosimmr_$output)$quantiles
        out_statistics[[i]] <- t(apply(out_all, 2, function(x) {
          return(c(mean = mean(x), sd = stats::sd(x)))
        }))
        # coda:::summary.mcmc.list(cosimmr_$output)$statistics[,1:2]
        out_cor[[i]] <- stats::cor(out_all)
        
        
        if ("quantiles" %in% type) {
          # Print out quantiles argument
          print(round(out_quantiles[[i]], 3))
        }
        
        if ("statistics" %in% type) {
          # Print out quantiles argument
          print(round(out_statistics[[i]], 3))
        }
        
        if ("correlations" %in% type) {
          # Print out quantiles argument
          print(round(out_cor[[i]], 3))
        }
        
        
        
        invisible(list(quantiles = out_quantiles, statistics = out_statistics, correlations = out_cor))
        
      }
    }
  } else {
    (return(message("incorrect object passed to function")))
  }
}


