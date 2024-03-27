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
#'  \code{\link{cosimmr_ffvb}}.
#' @param type The type of plot required. Can be one or more of 'isospace', 
#' 'betahistogram'
#' @param binwidth The width of the bins for the histogram. Defaults to 0.05
#' @param alpha The degree of transparency of the plots. Not relevant for
#' matrix plots
#' @param title The title of the plot.
#' @param individuals The individual number you wish to plot
#' @param covariates The covariate you wish to plot (for beta plots)
#' @param n_output The number of theta samples you wish to plot with. Defaults to 3600
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
plot.cosimmr_pred_out <-
  function(x,
           type = c("beta_histogram",
             "beta_boxplot",
             "p_ind"
           ),
           individuals = 1,
           covariates = 1,
           binwidth = 0.1,
           alpha = 0.5,
           title = NULL,
           n_output = 3600,
           ...) {
    if(inherits(x, "cosimmr_pred_out") == TRUE){
      title_input = title
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)
      
      # Iso-space plot is special as all groups go on one plot
      # Add in extra dots here as they can be sent to this plot function

      
      if("p_ind" %in% type){
        
        #Need to have a separate matrix for each ind value
        #So do all this in loop and repeat I think is easiest
        for(i in 1:length(individuals)){
          if(is.null(title_input) == TRUE){
            title = c(rep(NA, length(individuals)))
            for(j in 1:length(individuals)){
              title[j] = paste("p_ind plot: individual", individuals[j])
            }
          } else{title = rep(title_input, length(individuals))}
          curr_ind = individuals[i]
          out_all_p = x$p[curr_ind,,]
          
          
          colnames(out_all_p) = x$input$source_names
          
          df <- reshape2::melt(out_all_p)
          
          
          colnames(df) = c("Num", "Source", "Proportion")
          
          #add other plot types here maybe
          
          g <- ggplot(df, aes(
            x = Proportion
          )) +
            scale_fill_viridis(discrete = TRUE) +
            geom_histogram(binwidth = binwidth, alpha = alpha) +
            theme_bw() +
            ggtitle(title[i]) +
            facet_wrap("~ Source") +
            theme(legend.position = "none")
          print(g) 
          
          
        }
      }
      

      
      #Prep data
      #Data needs to be edited I think to make life easier
      for( i in 1:length(covariates)){
        
        
        
        beta = array(NA, dim = c(x$input$n_covariates, nrow(x$beta), x$input$n_sources))
        
        for(s in 1:nrow(x$theta)){
          for(k in 1:x$input$n_covariates){
            for(j in 1:x$input$n_sources){
              beta[k,s, j] = x$beta[s, (k-1)*x$input$n_sources + (j)]
            }
          }
        }
        
        #This is slightly wrong
        #I need to just think about it more
        
        out_all_beta = beta[i,,]
        colnames(out_all_beta) = x$input$source_names
        #I don't actually understand what this is doing
        df_beta <- reshape2::melt(out_all_beta)
        colnames(df_beta) = c("Num", "Source", "Beta")
        
        if("beta_histogram" %in% type){
          
          if(is.null(title_input) == TRUE){
            title = c(rep(NA, length(covariates)))
            for(j in 1:length(covariates)){
              title[j] = paste("beta histogram plot: covariate", covariates[j])
            }
          } else{title = rep(title_input, length(covariates))}
          
          #Histograms
          g <- ggplot(df_beta, aes(x = Beta)) +
            scale_fill_viridis(discrete = TRUE) +
            geom_histogram(binwidth = binwidth, alpha = alpha) +
            theme_bw() +
            ggtitle(title[i]) +
            facet_wrap("~ Source") +
            theme(legend.position = "none", 
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) +
            geom_vline(xintercept = 0, colour = "red")
          
          
          #Boxplot
          
          print(g)
        }
        
        
        
        
        if("beta_boxplot" %in% type){
          
          if(is.null(title_input) == TRUE){
            title = c(rep(NA, length(covariates)))
            for(j in 1:length(covariates)){
              title[j] = paste("beta boxplot: covariate", covariates[j])
            }
          } else{title = rep(title_input, length(covariates))}
          g <- ggplot(df_beta, aes(x = Beta)) +
            scale_fill_viridis(discrete = TRUE) +
            geom_boxplot() +
            theme_bw() +
            ggtitle(title[i]) +
            facet_wrap("~ Source")+
            geom_vline(xintercept = 0, colour = "red") +
            theme(axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
          
          print(g)
          
        }
        
        
        
      }
    }

    
    if (exists("g")) invisible(g) 
  }

