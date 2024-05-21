#' Plot the \code{cosimmr_input} data created from \code{cosimmr_load}
#'
#' This function creates iso-space (AKA tracer-space or delta-space) plots.
#' They are vital in determining whether the data are suitable for running in a
#' SIMM.
#'
#' It is desirable to have the vast majority of the mixture observations to be
#' inside the convex hull defined by the food sources. When there are more than
#' two tracers (as in one of the examples below) it is recommended to plot all
#' the different pairs of the food sources. See the vignette for further
#' details of richer plots.
#'
#' @param x An object created via the function \code{\link{cosimmr_load}}
#' @param tracers The choice of tracers to plot. If there are more than two
#' tracers, it is recommended to plot every pair of tracers to determine
#' whether the mixtures lie in the mixing polygon defined by the sources
#' @param title A title for the graph
#' @param xlab The x-axis label. By default this is assumed to be delta-13C but
#' can be made richer if required. See examples below.
#' @param ylab The y-axis label. By default this is assumed to be delta-15N in
#' per mil but can be changed as with the x-axis label
#' @param sigmas The number of standard deviations to plot on the source
#' values. Defaults to 1.
#' @param mix_name A optional string containing the name of the mixture
#' objects, e.g. Geese.
#' @param colour If TRUE (default) creates a plot. If not, puts the plot in
#' black and white
#' @param colour_by_cov if TRUE this allows users to colour the mixtures on the 
#' isospace plot by a specified covariate. Defaults to FALSE
#' @param  cov_name The name of the covariate the user wishes to colour the 
#' mixture points on the plot by
#' @param ggargs Extra arguments to be included in the ggplot (e.g. axis limits)
#' @param ...  Not used
#' 
#' @return isospace plot
#'
#' @import ggplot2
#' @import viridis
#' @import ggnewscale
#'
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
#' @seealso See \code{\link{plot.cosimmr_output}} for plotting the output of a
#' simmr run. See \code{\link{cosimmr_ffvb}} for running a cosimmr object once the
#' iso-space is deemed acceptable.
#' @examples
#' \dontrun{
#' # A simple example with 10 observations, 4 food sources and 2 tracers
#' data(geese_data_day1)
#' cosimmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ c(1,2,3,2,3,1,2,3,1),
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
#' ### A more complicated example with 30 obs, 3 tracers and 4 sources
#' data(simmr_data_2)
#' cosimmr_3 <- with(
#'   simmr_data_2,
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
#' # Plot 3 times - first default d13C vs d15N
#' plot(cosimmr_3)
#' # Now plot d15N vs d34S
#' plot(cosimmr_3, tracers = c(2, 3))
#' # and finally d13C vs d34S
#' plot(cosimmr_3, tracers = c(1, 3))
#' }
#' @export
plot.cosimmr_input <-
  function(x,
           tracers = c(1, 2),
           title = "Tracers plot",
           xlab = colnames(x$mixtures)[tracers[1]],
           ylab = colnames(x$mixtures)[tracers[2]],
           sigmas = 1,
           mix_name = "Mixtures",
           colour = TRUE,
           colour_by_cov = FALSE,
           cov_name = NULL,
           ggargs = NULL,
           ...) {
    
    #This selects the correct column from the covariates df to colour by
    if(colour_by_cov == TRUE){
      cov_selected_col = subset(as.matrix(x$covariates_df), select =  cov_name)
    }
    
    curr_mix <- x$mixtures 
    
    
    # First get the mean corrected sources and the sd corrected sources
    source_means_c <- x$source_means + x$correction_means
    source_sds_c <- sqrt(x$source_sds^2 + x$correction_sds^2)
    
    # Set up data frame for ggplot - have to do it this stupid way because of cran
    x2 <- unlist(c(source_means_c[, tracers[1]], curr_mix[, tracers[1]]))
    x_lower <- unlist(c(source_means_c[, tracers[1]] - sigmas * source_sds_c[, tracers[1]], curr_mix[, tracers[1]]))
    x_upper <- unlist(c(source_means_c[, tracers[1]] + sigmas * source_sds_c[, tracers[1]], curr_mix[, tracers[1]]))
    
    if (ncol(curr_mix) > 1) {
      y <- unlist(c(source_means_c[, tracers[2]], curr_mix[, tracers[2]]))
      y_lower <- unlist(c(source_means_c[, tracers[2]] - sigmas * source_sds_c[, tracers[2]], curr_mix[, tracers[2]]))
      y_upper <- unlist(c(source_means_c[, tracers[2]] + sigmas * source_sds_c[, tracers[2]], curr_mix[, tracers[2]]))
    }
    
    
    #Set the sources column - okay for factors but need something else for numeric
    if(colour_by_cov == FALSE){
      Source <- factor(c(x$source_names, rep(mix_name, nrow(curr_mix))), levels = c(mix_name, x$source_names))
    } else {
      if(is.numeric(as.matrix(cov_selected_col))){
        Source <- c(x$source_names, (as.matrix(cov_selected_col)))#factor(c(x$source_names, rep(mix_name, nrow(curr_mix))), levels = c(mix_name, x$source_names))
      } else{
        Source <- factor(c(x$source_names, (as.matrix(cov_selected_col))))
      }
    }
    
    
    size <- c(rep(0.5, x$n_sources), rep(0.5, nrow(curr_mix)))
    
    if (ncol(curr_mix) == 1) {
      df <- data.frame(x = x2, x_lower, x_upper, Source, size, y = Source)
    } else {
      df <- data.frame(x = x2, y = y, x_lower, y_lower, x_upper, y_upper, Source, size)
    }
    
    
    # Plot for bivariate mixtures
    if (ncol(curr_mix) > 1) {
      if (colour) {
        if(colour_by_cov == FALSE){
          g <- ggplot(data = df, aes(x = x, y = y, colour = Source)) +
            scale_color_viridis(discrete = TRUE) +
            theme_bw() +
            labs(x = xlab, y = ylab, title = title) +
            geom_errorbarh(aes(xmax = x_upper, xmin = x_lower, height = 0)) +
            geom_pointrange(aes(x = x, y = y, ymax = y_upper, ymin = y_lower, shape = Source)) +
            scale_shape_manual(values = 1:nlevels(df$Source)) +
            theme(legend.title = element_blank(), legend.key = element_blank()) +
            guides(color = guide_legend(override.aes = list(linetype = c(rep(0, 1), rep(1,(nlevels(df$Source)-1)))))) +
            ggargs
        }
        else if(is.numeric(cov_selected_col)){
          
          # split data$Source into mixtures and sources
          group_data <- df[-c(1:x$n_sources), ]
          non_group_data <- df[c(1:x$n_sources), ]
          
          # Generate Shapes ---------------------------------------------------------
          
          # This function sets the 'Groups' to be circles and anything else as different shapes.
          # This might be a bit overkill, but i thought it would be better to automate it rather than hardcoding it
          generate_shapes <- function(data_Source) {
            # initial shapes
            non_group_shapes <- 1:(x$n_sources) # not sure how many you'll need here or if there will always be 3
            
            # set up vector to store shape
            shapes_vector <- integer(length(data_Source))
            
            # give 'group' a circle shape
            is_group <- c(rep(FALSE, x$n_sources), rep(TRUE, x$n_obs))
            shapes_vector[is_group] <- x$n_obs
            
            # give non-group other shapes (wrap around if there are more non-groups than selected shapes)
            non_group_indices <- which(!is_group)
            shapes_vector[non_group_indices] <- non_group_shapes[(non_group_indices - 1) %% length(non_group_shapes) + 1]
            names(shapes_vector) <- data_Source
            
            return(shapes_vector)
          }
          
          # get shapes for the plot
          shapes <- generate_shapes(data_Source = df$Source)
          
          # create dummy data for 'Mixtures'. This is only used to include 'Mixture' in the legend
          non_group_data_extended <- rbind(non_group_data, data.frame(
            x = NA,
            y = NA,
            y_upper = NA,
            y_lower = NA,
            x_upper = NA,
            x_lower = NA,
            Source = "Mixtures",
            size = NA
          ))
          
          
          #cols vector
          
          
          #Going to do this in a v hacky way
          al_order = sort(non_group_data_extended$Source) #Sorting in alphabetical order to match ggplot legend
          int_mix = grep("Mixtures", al_order, value = FALSE) #Finding mixtures
          
          cols = c(rep(NA, x$n_sources + 1)) #now creating a colours column
          cols = viridis::mako(x$n_sources + 1) #Filling it with viridis pallette
          
          cols[int_mix] = c("Mixtures" = "black") #Making the right colour black for mixtures circle
          group_data$Source = as.numeric(group_data$Source) #Was getting an error otherwise - somewhere along the line this converts to chr
          
          g <- ggplot() +
            theme_bw() +
            geom_errorbarh(data = non_group_data, aes(y = y, xmax = x_upper, xmin = x_lower, height = 0, color = Source)) +
            geom_pointrange(data = non_group_data, aes(x = x, y = y, ymax = y_upper, ymin = y_lower, shape = Source, color = Source)) +
            geom_point(
              data = non_group_data_extended[non_group_data_extended$Source == "Mixtures", ],
              aes(x = x, y = y, shape = Source, color = Source), size = 3
            ) + 
       #     theme(legend.title = element_blank(), legend.key = element_blank())+
            scale_color_manual(  # Mixtures point
              name = "",
              values = cols
            ) +
            scale_shape_manual(
              name = "",
              values = c(shapes, "Mixtures" = 16)
            ) +
            ggnewscale::new_scale_color() +
            geom_point(data = group_data, aes(x = x, y = y, color = Source), size = 3) +
            scale_color_viridis_c(
              name = cov_name,
              guide = guide_colorbar(
                order = 2,
                frame.colour = "black",
                ticks.colour = "black"
              )
            ) + labs(x = xlab, y = ylab, title = title) 
          
          # due to the way im adding Mixtures as a kind of non-existent point, this will produce a warning saying:
          # Removed 1 row containing missing values or values outside the scale range. To suppress ggplot warnings, you gotta 
          # wrap them in a print()
        #  suppressWarnings(print(p))
          
          
        }else if(!is.numeric(cov_selected_col)){
          #Need to pull out the groups vs sources I think? To remove the lines on the sides
          group_data <- df[-c(1:x$n_sources), ]
          non_group_data <- df[c(1:x$n_sources), ]
          
         
          
          g <- ggplot() +
            theme_bw() +
            geom_errorbarh(data = non_group_data, aes(y = y, xmax = x_upper, xmin = x_lower, height = 0, color = Source)) +
            geom_pointrange(data = non_group_data, aes(x = x, y = y, ymax = y_upper, ymin = y_lower, shape = Source, color = Source)) +
            geom_point(data = group_data, aes(x = x, y = y, shape = Source, color = Source), size = 3) +
            scale_colour_viridis(discrete = TRUE) +  labs(x = xlab, y = ylab, title = title) + 
            scale_shape_manual(values = c(1:nlevels(df$Source))) +
            theme(legend.title = element_blank(), legend.key = element_blank())
          

            
        }
          
        
      } else {
        g <- ggplot(data = df, aes(x = x, y = y, colour = Source)) +
          theme_bw() +
          labs(x = xlab, y = ylab, title = title) +
          geom_errorbarh(aes(xmax = x_upper, xmin = x_lower, height = 0)) +
          # geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
          geom_pointrange(aes(x = x, y = y, ymax = y_upper, ymin = y_lower, shape = Source)) +
          scale_shape_manual(values = 1:nlevels(df$Source)) +
          theme(legend.title = element_blank(), legend.key = element_blank()) +
          guides(color = guide_legend(override.aes = list(linetype = c(rep(0, 1), rep(1, x$n_sources))))) +
          scale_colour_grey() +
          ggargs  +
          theme(legend.title = element_blank(), legend.key = element_blank())
      }
    }
    
    # Plot for univariate mixtures
    if (ncol(curr_mix) == 1) {
      if (colour) {
        g <- ggplot(data = df, aes(x = x, y = y, colour = Source)) +
          scale_color_viridis(discrete = TRUE) +
          theme_bw() +
          theme(axis.title.y = element_blank()) +
          labs(x = xlab, title = title) +
          geom_errorbarh(aes(xmax = x_upper, xmin = x_lower, height = 0)) +
          # geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
          geom_point(aes(shape = Source)) +
          scale_shape_manual(values = 1:nlevels(df$Source)) +
          theme(legend.position = "None") +
          guides(color = guide_legend(override.aes = list(linetype = c(rep(0, 1), rep(1, x$n_sources))))) +
          ggargs  +
          theme(legend.title = element_blank(), legend.key = element_blank())
      } else {
        g <- ggplot(data = df, aes(x = x, y = y, colour = Source)) +
          scale_color_grey() +
          theme_bw() +
          theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          ) +
          labs(x = xlab, title = title) +
          geom_errorbarh(aes(xmax = x_upper, xmin = x_lower, height = 0)) +
          # geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
          geom_point(aes(shape = Source)) +
          scale_shape_manual(values = 1:nlevels(df$Source)) +
          theme(legend.title = element_blank(), legend.key = element_blank()) +
          guides(color = guide_legend(override.aes = list(linetype = c(rep(0, 1), rep(1, x$n_sources))))) +
          ggargs  +
          theme(legend.title = element_blank(), legend.key = element_blank())
      }
    }
   #  if(colour_by_cov == TRUE){
   #      if(is.numeric(cov_selected_col)){
   #    suppressWarnings(print(g))
   #    invisible(g)
   #  }
   #  }
   #  if(colour_by_cov == FALSE){ 
   #  print(g)
   #  invisible(g)
   # }
    suppressWarnings(print(g))
    invisible(g)
  }
