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
             "beta_histogram",
             "beta_boxplot",
             "p_ind",
             "p_mean"
             ),
           ind = 1,
           covariates = 1,
           binwidth = 0.1,
           alpha = 0.5,
           title = NULL,
           n_output = 3600,
           ...) {
    if(inherits(x, "cosimmr_output") == TRUE){
      title_input = title
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)
      
      # Iso-space plot is special as all groups go on one plot
      # Add in extra dots here as they can be sent to this plot function
      if ("isospace" %in% type) {
        if(is.null(title_input) == TRUE){
          title = "isospace plot"
        } else{title = title_input}
        graphics::plot(x$input, title = title, ...)
        
      }
      
      if("p_ind" %in% type){
        
        #Need to have a separate matrix for each ind value
        #So do all this in loop and repeat I think is easiest
        for(i in 1:length(ind)){
          if(is.null(title_input) == TRUE){
            title = c(rep(NA, length(ind)))
            for(j in 1:length(ind)){
              title[j] = paste("p_ind plot: individual", ind[j])
            }
          } else{title = rep(title_input, length(ind))}
          curr_ind = ind[i]
        out_all_p = x$output$BUGSoutput$sims.list$p[curr_ind,,]
        

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
      
      if("p_mean" %in% type){
        #I think I'm right in saying that the mean should be
        #When all the x's equal zero because its mean centred
        #Also need to check for intercept I guess?
        
        if(is.null(title_input) == TRUE){
            title = "p_mean_plot"
          } else{title = title_input}
       

      #So check if intercept = TRUE and then print a vector
        if(x$input$intercept == TRUE){
          x_pred = c(1, rep(0, (ncol(x$input$x_scaled) - 1)))
        } else if(x$input$intercept == FALSE){
          x_pred = c(rep(0, (ncol(x$input$x_scaled))))
        }

        thetares= x$output$theta
        K = x$input$n_sources
        n_tracers = x$input$n_tracers
        n_covariates = ncol(x$input$x_scaled)



        sigma <- (1/sqrt(thetares[,(K*n_covariates + 1):(K*n_covariates + n_tracers)]))

        #p_sample = array(NA, dim = c(1, n_output, K))
        p_sample = matrix(ncol = K, nrow = n_output)

        beta = array(thetares[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))

        f <- array(NA, dim = c(1, K, n_output))

        for(s in 1:n_output){
          f[,,s] = (x_pred) %*% beta[s,,]
        }

        for(j in 1:n_output){
           # p_sample[1,j, ] 
         p_sample[j,] <- exp(f[1,1:K, j]) / (sum((exp(f[1,1:K, j]))))
        }
        
        colnames(p_sample) = x$input$source_names

        df_p_mean <- reshape2::melt(p_sample)


        colnames(df_p_mean) = c("Num", "Source", "Proportion")

        g <- ggplot(df_p_mean, aes(
          x = Proportion
        )) +
          scale_fill_viridis(discrete = TRUE) +
          geom_histogram(binwidth = binwidth, alpha = alpha) +
          theme_bw() +
          ggtitle(title) +
          facet_wrap("~ Source") +
          theme(legend.position = "none")
        print(g)



         }
  
      #Prep data
      #Data needs to be edited I think to make life easier
      
    for( i in 1:length(covariates)){
      out_all_beta = x$output$beta[,i,]
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

