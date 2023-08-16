#' Plot different features of an object created from  \code{\link{cosimmr_ffvb}}.
#'
#' This function allows for 4 different types of plots of the simmr output
#' created from \code{\link{cosimmr_ffvb}}. The 
#' types are: plot of beta values
#'
#' The matrix plot should form a necessary part of any SIMM analysis since it
#' allows the user to judge which sources are identifiable by the model.
#' Further detail about these plots is provided in the vignette.
#'
#' @param x An object of class \code{cosimmr_output} created via
#'  \code{\link{simmr_ffvb}}.
#' @param type The type of plot required. Can be one or more of 'isospace', 
#' 'betahistogram'
#' @param binwidth The width of the bins for the histogram. Defaults to 0.05
#' @param alpha The degree of transparency of the plots. Not relevant for
#' matrix plots
#' @param title The title of the plot.
#' @param ggargs Extra arguments to be included in the ggplot (e.g. axis limits)
#' @param ...  Currently not used
#'
#' @import ggplot2
#' @import graphics
#' @import viridis
#' @importFrom reshape2 "melt"
#' @importFrom stats "cor"
#'
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See  \code{\link{cosimmr_ffvb}} for 
#' creating objects suitable for this function, and many more examples. See 
#' also \code{\link{cosimmr_load}} for creating simmr objects, 
#' \code{\link{plot.cosimmr_input}} for creating isospace plots.
#' 
#' @examples
#'
#' \dontrun{
#' # A simple example with 10 observations, 2 tracers and 4 sources
#'
#' # The data
#' data(geese_data_day1)
#'
#' # Load into simmr
#' simmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ c(1,2,3,3,2,3,2,1,2),
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means
#'   )
#' )
#' # Plot
#' plot(simmr_1)
#'
#'
#' # FFVB run
#' simmr_1_out <- cosimmr_ffvb(simmr_1)
#'
#'plot(simmr_1_out, type = c("isospace", "beta_hist"))
#' }
#' @export
plot.cosimmr_output <-
  function(x,
           type = c(
             "isospace",
             "beta_hist"
             ),
           covariates = 1,
           binwidth = 10,
           alpha = 0.5,
           title = if (length(covariates) == 1) {
             "simmr output plot"
           } else {
             paste("simmr output plot: covariate", covariates)
           },
           ggargs = NULL,
           ...) {
    if(inherits(x, "cosimmr_output") == TRUE){
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)
      
      # Iso-space plot is special as all groups go on one plot
      # Add in extra dots here as they can be sent to this plot function
      if ("isospace" %in% type) {
        graphics::plot(x$input, title = title, ...)
      }
      
      #Prep data
      #Data needs to be edited I think to make life easier
      
    for( i in 1:length(covariates)){
      out_all_beta = x$output$beta[,i,]
      #I don't actually understand what this is doing
      df <- reshape2::melt(out_all_beta)
      colnames(df) = c("Num", "Source", "Beta")
      
      if("beta_hist" %in% type){
        
        #Histograms
        g <- ggplot(df, aes(
          x = Beta
        )) +
          scale_fill_viridis(discrete = TRUE) +
          geom_histogram(binwidth = binwidth, alpha = alpha) +
          theme_bw() +
          ggtitle(title[i]) +
          facet_wrap("~ Source") +
          theme(legend.position = "none") +
          ggargs +
          geom_vline(xintercept = 0, colour = "red")
        
        
        #Boxplot
        g <- ggplot(df, aes(x = Beta)) +
          scale_fill_viridis(discrete = TRUE) +
          geom_boxplot() +
         theme_bw() +
        ggtitle(title[i]) +
         facet_wrap("~ Source")+
  geom_vline(xintercept = 0, colour = "red") +
          ggargs

        print(g)
      }
      
      
      if("x" %in% type){
        #add other plot types here maybe
      }
      
      
      
    }
    }
    if (exists("g")) invisible(g) 
  }
