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
#' 'beta_histogram', 'beta_boxplot', 'prob_histogram', 'prob_density', 'covariates_plot'
#' @param binwidth The width of the bins for the histogram. Defaults to 0.05
#' @param alpha The degree of transparency of the plots. Not relevant for
#' matrix plots
#' @param obs The observation number you wish to plot
#' @param cov_name The name of the covariate you wish to plot (for beta and covariates plot)
#' @param title The title of the plot.
#' @param n_output The number of theta samples you wish to plot with. Defaults to 3600
#' @param source The number or name of the source you wish to plot over for 
#' 'covariates_plot', defaults to NULL which means all sources are used
#' @param one_plot Whether to plot line covariates plot on one plot. Defaults to FALSE
#' @param ...  Currently not used
#'
#' @import ggplot2
#' @import graphics
#' @import viridis
#' @importFrom reshape2 "melt"
#' @importFrom stats "cor"
#'
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
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
#'     formula = mixtures ~ 1,
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
             "prop_histogram",
             "prop_density",
             "covariates_plot"
             ),
           obs = 1,
           cov_name = NULL, 
           binwidth = 0.05,
           alpha = 0.5,
           title = NULL,
           n_output = 3600,
           source = NULL,
           one_plot = FALSE,
           ...) {
    if(inherits(x, "cosimmr_output") == TRUE){
      title_input = title
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)
      
      #want to extract the right covariate number using the name?
      n_cov = length(cov_name)
      covariates = c(rep(NA, n_cov))
      if(ncol(x$input$original_x) ==1 & x$input$intercept == TRUE  & (("beta_histogram" %in% type)|("beta_boxplot" %in% type)|("covariates_plot" %in% type))){
        message("You have selected a plot type that requires covariates but 
the model does not contain covariates. Please reselect
model or plots to create and rerun.")
      } else{
      
      
      
      # Iso-space plot is special as all groups go on one plot
      # Add in extra dots here as they can be sent to this plot function
      if ("isospace" %in% type) {
        if(is.null(title_input) == TRUE){
          title = "isospace plot"
        } else{title = title_input}
        graphics::plot(x$input, title = title, ...)
        
      }
      
      
      
    
        
        #Need to have a separate matrix for each ind value
        #So do all this in loop and repeat I think is easiest
        for(i in 1:length(obs)){
          if(is.null(title_input) == TRUE){
            title = c(rep(NA, length(obs)))
            for(j in 1:length(obs)){
              title[j] = paste("Proportions: Observation", obs[j])
            }
          } else{title = rep(title_input, length(obs))}
          curr_ind = obs[i]
        out_all_p = x$output$BUGSoutput$sims.list$p[curr_ind,,]
        

         colnames(out_all_p) = x$input$source_names

        df <- reshape2::melt(out_all_p)
        
        
        colnames(df) = c("Num", "Source", "Proportion")
        
        #add other plot types here maybe
        if("prop_histogram" %in% type){
          g <- ggplot(df, aes(
            x = Proportion,
            fill = Source
          )) +
            scale_fill_viridis(discrete = TRUE) +
            geom_histogram(binwidth = binwidth, alpha = alpha) +
            theme_bw() +
            ggtitle(title[i]) +
            facet_wrap("~ Source") +
            theme(legend.position = "none")
          print(g) 
          
          
        }
      
      
      if("prop_density" %in% type){
        g <- ggplot(df, aes(
          x = Proportion,
          fill = Source
        )) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(aes(y = after_stat(density)), alpha = alpha, linetype = 0) +
          theme_bw() +
          theme(legend.position = "none") +
          ggtitle(title[i]) +
          ylab("Density") +
          facet_wrap("~ Source")
        print(g)
      }
        }
      
      
  
      #Prep data
      #Data needs to be edited I think to make life easier
      
    for(l in 1:length(covariates)){
      if("beta_histogram" %in% type){
   
      if(x$input$intercept == FALSE){
        for(i in 1:n_cov){
          covariates[i] = grep(paste0(cov_name[i]), colnames(x$input$covariates_df), value = FALSE)
        }
      }else if(x$input$intercept == TRUE){
        for(i in 1:n_cov){
          covariates[i] = grep(paste0(cov_name[i]), colnames(x$input$covariates_df), value = FALSE)
        }
      }
      
      
      if(x$input$intercept == TRUE){
      cov_ind =covariates[l] +1} else{
        cov_ind =covariates[l]
      }

      beta = array(NA, dim = c(x$input$n_covariates, nrow(x$output$beta), x$input$n_sources))

      for(s in 1:nrow(x$output$theta)){
      for(c in 1:x$input$n_covariates){
        for(j in 1:x$input$n_sources){
          beta[c,s, j] = x$output$beta[s, (c-1)*x$input$n_sources + (j)]
        }
      }
      }



      out_all_beta = beta[cov_ind,,]
      colnames(out_all_beta) = x$input$source_names
      #I don't actually understand what this is doing
      df_beta <- reshape2::melt(out_all_beta)
      colnames(df_beta) = c("Num", "Source", "Beta")

     
 
          if(is.null(title_input) == TRUE){
            title = c(rep(NA, length(covariates)))
            for(c in 1:length(covariates)){
              title[c] = paste("beta histogram plot: covariate", cov_name[c])
            }
          } else{title = rep(title_input, length(covariates))}
        

        #Histograms
        g <- ggplot(df_beta, aes(x = Beta)) +
          scale_fill_viridis(discrete = TRUE) +
          geom_histogram(binwidth = binwidth, alpha = alpha) +
          theme_bw() +
          ggtitle(title[l]) +
          facet_wrap("~ Source") +
          theme(legend.position = "none",
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
          geom_vline(xintercept = 0, colour = "red")


        #Boxplot

        print(g)
      } 

      if("beta_boxplot" %in% type){
        
        if(x$input$intercept == FALSE){
          for(i in 1:n_cov){
            covariates[i] = grep(paste0(cov_name[i]), colnames(x$input$covariates_df), value = FALSE)
          }
        }else if(x$input$intercept == TRUE){
          for(i in 1:n_cov){
            covariates[i] = grep(paste0(cov_name[i]), colnames(x$input$covariates_df), value = FALSE)
          }
        }
        
        
        if(x$input$intercept == TRUE){
          cov_ind =covariates[l] +1} else{
            cov_ind =covariates[l]
          }
        
        beta = array(NA, dim = c(x$input$n_covariates, nrow(x$output$beta), x$input$n_sources))
        
        for(s in 1:nrow(x$output$theta)){
          for(c in 1:x$input$n_covariates){
            for(j in 1:x$input$n_sources){
              beta[c,s, j] = x$output$beta[s, (c-1)*x$input$n_sources + (j)]
            }
          }
        }
        
        
        
        out_all_beta = beta[cov_ind,,]
        colnames(out_all_beta) = x$input$source_names
        #I don't actually understand what this is doing
        df_beta <- reshape2::melt(out_all_beta)
        colnames(df_beta) = c("Num", "Source", "Beta")
        

        if(is.null(title_input) == TRUE){
          title = c(rep(NA, length(covariates)))
          for(c in 1:length(covariates)){
            title[c] = paste("beta boxplot: covariate", cov_name[c])
          }
        } else{title = rep(title_input, length(covariates))}
        g <- ggplot(df_beta, aes(x = Beta)) +
          scale_fill_viridis(discrete = TRUE) +
          geom_boxplot() +
          theme_bw() +
          ggtitle(title[l]) +
          facet_wrap("~ Source")+
          geom_vline(xintercept = 0, colour = "red") +
          theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())

        print(g)

      }
      
      if("covariates_plot" %in% type){
        if(x$input$intercept == FALSE){
          for(i in 1:n_cov){
            covariates[i] = grep(paste0(cov_name[i]), colnames(x$input$covariates_df), value = FALSE)
          }
        }else if(x$input$intercept == TRUE){
          for(i in 1:n_cov){
            covariates[i] = grep(paste0(cov_name[i]), colnames(x$input$covariates_df), value = FALSE)
          }
        }
        
        
        if(x$input$intercept == TRUE){
          cov_ind =covariates[l] +1} else{
            cov_ind =covariates[l]
          }
        

        #want to let them choose the food source by name?? I guess??
        if(is.null(source)){
          source_loop = x$input$n_sources
          source_chosen = c(x$input$source_names)
        }else if(length(source) == 1){
          source_loop = 1
          source_chosen = source
        }else{
          source_loop = length(source)
          source_chosen = source
          }
        
        if(one_plot == TRUE){
          
          #so now we have the source name (s_name) and the number of the source (source_n)
          n_samples = length(x$output$BUGSoutput$sims.list$p[1,,1])
          n_ind = x$input$n_obs
          
          #check which plot we want to make - boxplot or lineplot
          box_or_line = NULL
          
          
          cov_selected_col = matrix(c((x$input$covariates_df)[,l]), ncol = 1)
          
          colnames(cov_selected_col) = colnames(x$input$covariates_df)[l]
          
          #Just chooses which plot we're making
          if(is.numeric(cov_selected_col)){
            box_or_line = "LINE"
          }else {
            box_or_line = "BOX"
          }
          
          if(box_or_line == "LINE"){
            n_sources = x$input$n_sources
            line_array = array(NA, dim = c(n_ind, n_samples, n_sources))
            mean_line_mat = matrix(NA, nrow = n_ind, ncol = n_sources)
            save_sd =  matrix(NA, nrow = n_ind, ncol = n_sources)
            for(s in 1:n_sources){
            line_array[,,s] = matrix(x$output$BUGSoutput$sims.list$p[,,s], nrow = n_ind, ncol = n_samples)
            mean_line_mat[,s] = rowMeans(line_array[,,s])
            for(i in 1:n_ind){
              save_sd[i,s] = sd(line_array[i,,s])
            }
            }
            
            out_all_mean <- mean_line_mat
            colnames(out_all_mean) <- x$input$source_names
            df_mean <- reshape2::melt(out_all_mean)
            
            colnames(df_mean) <- c("Num", "Source", "Mean")
            
            out_all_sd <- save_sd
            colnames(out_all_sd) <- x$input$source_names
            df_sd <- reshape2::melt(out_all_sd)
            
            colnames(df_sd) <- c("Num", "Source", "SD")
            
            df_plot = data.frame(mean = df_mean$Mean,
                                 sd = df_sd$SD,
                                 cov = x$input$covariates_df[,l],
                                 psd = df_mean$Mean + df_sd$SD,
                                 nsd = df_mean$Mean - df_sd$SD,
                                 num = df_mean$Num,
                                 Source = df_mean$Source)
            

            
            g <- ggplot(data = df_plot, aes(x = cov, y = mean, colour = Source)) +
              geom_ribbon(data = df_plot, aes(ymin = nsd, ymax = psd, fill = Source), alpha = alpha) +
              geom_line() + xlab(colnames(cov_selected_col)) +
              ylab("Proportion (± 1sd)") + ggtitle(paste0("Change in consumption over ", colnames(cov_selected_col)))
            
            print(g)
            
            
          }else if(box_or_line == "BOX"){
            message("Cannot facet wrap boxplot")
          
          }
        }
        else if(one_plot == FALSE){
        for(s in 1:source_loop){
         
        
         if(is.numeric(source)){
          source_n = source_chosen[s]
          s_name = x$input$source_names[source_n]
          
        }else{
          source_n = grep(source_chosen[s], x$input$source_names, value = FALSE)
          s_name = source_chosen[s]
        }
          
#so now we have the source name (s_name) and the number of the source (source_n)
        n_samples = length(x$output$BUGSoutput$sims.list$p[1,,1])
        n_ind = x$input$n_obs
        
        #check which plot we want to make - boxplot or lineplot
        box_or_line = NULL


          cov_selected_col = matrix(c((x$input$covariates_df)[,l]), ncol = 1)

          colnames(cov_selected_col) = colnames(x$input$covariates_df)[l]
        
        #Just chooses which plot we're making
        if(is.numeric(cov_selected_col)){
          box_or_line = "LINE"
        }else {
          box_or_line = "BOX"
        }
        
        if(box_or_line == "LINE"){
          line_mat = matrix(x$output$BUGSoutput$sims.list$p[,,s], nrow = n_ind, ncol = n_samples)
          mean_line_mat = rowMeans(line_mat)
          save_sd = c(rep(NA, n_ind))
          for(i in 1:n_ind){
            save_sd[i] = sd(line_mat[i,])
          }
          
          df_plot = data.frame(mean = mean_line_mat,
                               sd = save_sd,
                               cov = x$input$covariates_df[,l],
                              psd = (mean_line_mat + save_sd),
                              nsd = (mean_line_mat - save_sd))
          

          
          
          g <- ggplot(data = df_plot, aes(x = cov, y = mean)) +
            geom_ribbon(data = df_plot, aes(ymin = nsd, ymax = psd), alpha = alpha) +
            geom_line() + 
            ggtitle(paste0("Proportion changing for ", s_name, " consumption over ", colnames(cov_selected_col))) +
            xlab(paste0(colnames(cov_selected_col))) + ylab("Proportion")
            
            print(g)
          
          
        }else if(box_or_line == "BOX"){
          #First make this the length of each group
          n_groups = length(unique(cov_selected_col))
          g_names = c(unique(cov_selected_col))
          #group_vec = c(rep(NA, n_groups))
          grep_values = c(rep(NA, n_groups))
          
          for(i in 1:n_groups){
            grep_values[i] = grep(g_names[i], cov_selected_col)[1]
          }
          
          
          #Now we want to extract an individual for each group
          
          
          # l is the covariate number
         # a = matrix(x$output$BUGSoutput$sims.list$p[grep_values,,l], ncol = n_samples)
          ind_mat = matrix(nrow = n_samples, ncol = n_groups)
          
          ########MISTAKE HERE!!!!!!!!!!!!!!!!!!!!!!
           for(i in 1:n_groups){
             ind_val = grep_values[i]
             ind_mat[,i] = x$output$BUGSoutput$sims.list$p[ind_val,,s]
           }
          
          ind_vec =c(ind_mat)
          
          #Now we just have to make the group names repeat n_samples times each
          g_names_mat = matrix(nrow = n_samples, ncol = n_groups)
          for(i in 1:n_groups){
            g_names_mat[,i] = c(rep(paste0(g_names[i]), n_samples))
          }
          
          g_names_rep_vec = c(g_names_mat)
          
          df_plot = data.frame(samples = ind_vec, Group = g_names_rep_vec)
          
          
          g = ggplot(data = df_plot, aes(x = Group, y = samples, colour = Group))  + 
            geom_boxplot() + 
            ggtitle(paste0(s_name, " consumption over ",colnames(cov_selected_col), " covariate")) +
            xlab(paste0(colnames(cov_selected_col))) + ylab("Proportion")
          
          print(g)
          
        }
     
        }
      
      
      
      
    }
    
      } 
      
      
      
      
      
      } #cov loop bracket
    
    if (exists("g")) invisible(g) 
      }
    }
    
  }

