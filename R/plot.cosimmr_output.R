#This needs to be edited - not sure how we want to go about plots


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
#' @param group Which groups to plot. Can be a single group or multiple groups
#' @param mix_name A optional string containing the name of the mixture
#' objects, e.g. Geese.
#' @param colour If TRUE (default) creates a plot. If not, puts the plot in
#' black and white
#' @param ggargs Extra arguments to be included in the ggplot (e.g. axis limits)
#' @param ...  Not used
#' 
#' @return isospace plot
#'
#' @import ggplot2
#' @import viridis
#'
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{plot.simmr_output}} for plotting the output of a
#' simmr run. See \code{\link{simmr_mcmc}} for running a simmr object once the
#' iso-space is deemed acceptable.
#' @examples
#'
#' # A simple example with 10 observations, 4 food sources and 2 tracers
#' data(geese_data_day1)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   simmr_load(
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
#' ### A more complicated example with 30 obs, 3 tracers and 4 sources
#' data(simmr_data_2)
#' simmr_3 <- with(
#'   simmr_data_2,
#'   simmr_load(
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
#' # Plot 3 times - first default d13C vs d15N
#' plot(simmr_3)
#' # Now plot d15N vs d34S
#' plot(simmr_3, tracers = c(2, 3))
#' # and finally d13C vs d34S
#' plot(simmr_3, tracers = c(1, 3))
#' # See vignette('simmr') for fancier x-axis labels
#'
#' # An example with multiple groups - the Geese data from Inger et al 2006
#' data(geese_data)
#' simmr_4 <- with(
#'   geese_data,
#'   simmr_load(
#'     mixtures = mixtures,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means,
#'     group = groups
#'   )
#' )
#'
#' # Print
#' simmr_4
#'
#' # Plot
#' plot(simmr_4,
#'   xlab = expression(paste(delta^13, "C (\u2030)", sep = "")),
#'   ylab = expression(paste(delta^15, "N (\u2030)", sep = "")),
#'   title = "Isospace plot of Inger et al Geese data"
#' ) #'
#' @export
plot.cosimmr_output <-
  function(x,
           tracers = c(1, 2),
           title = "Tracers plot",
           xlab = colnames(x$input$mixtures)[tracers[1]],
           ylab = colnames(x$inputmixtures)[tracers[2]],
           sigmas = 1,
           group = 1:x$input$n_groups,
           mix_name = "Mixtures",
           ggargs = NULL,
           colour = TRUE,
           ...) {
    
    if (inherits(x, "cosimmr_output") == TRUE) {
    
    # Get mixtures to match current group(s)
    curr_rows <- which(x$input$group_int %in% group)
    curr_mix <- x$input$mixtures[curr_rows, , drop = FALSE]
    curr_n_groups <- length(group)
    
    # First get the mean corrected sources and the sd corrected sources
    source_means_c <- x$input$source_means + x$input$correction_means
    source_sds_c <- sqrt(x$input$source_sds^2 + x$input$correction_sds^2)
    
    # Set up data frame for ggplot - have to do it this stupid way because of cran
    x2 <- unlist(c(source_means_c[, tracers[1]], curr_mix[, tracers[1]]))
    x_lower <- unlist(c(source_means_c[, tracers[1]] - sigmas * source_sds_c[, tracers[1]], curr_mix[, tracers[1]]))
    x_upper <- unlist(c(source_means_c[, tracers[1]] + sigmas * source_sds_c[, tracers[1]], curr_mix[, tracers[1]]))
    
    if (ncol(curr_mix) > 1) {
      y <- unlist(c(source_means_c[, tracers[2]], curr_mix[, tracers[2]]))
      y_lower <- unlist(c(source_means_c[, tracers[2]] - sigmas * source_sds_c[, tracers[2]], curr_mix[, tracers[2]]))
      y_upper <- unlist(c(source_means_c[, tracers[2]] + sigmas * source_sds_c[, tracers[2]], curr_mix[, tracers[2]]))
    }
    
    if (x$input$n_groups == 1) {
      Source <- factor(c(x$input$source_names, rep(mix_name, nrow(curr_mix))), levels = c(mix_name, x$input$source_names))
      p_means = matrix(NA, nrow = (x$input$n_obs + x$input$n_sources), ncol = x$input$n_sources)
      # for(k in 1:x$input$n_sources){
      #   p_means[k,] = c(rep(1, x$input$n_sources))
      # }
      p_means[1:x$input$n_sources,] = diag(x$input$n_sources)
      for(k in 1:x$input$n_obs){
        p_means[k+x$input$n_sources,] = c(colMeans(x$output[[group[i]]]$BUGSoutput$sims.list$p[k,,]))
      }
      colnames(p_means) = x$input$source_names
    } else {
      Source <- factor(c(x$input$source_names, paste(mix_name, x$input$group[curr_rows])), levels = c(paste(mix_name, unique(x$input$group[curr_rows])), x$source_names))
      p_means = matrix(NA, nrow = (x$input$n_obs + x$input$n_sources), ncol = x$input$n_sources)
      # for(k in 1:x$input$n_sources){
      #   p_means[k,] = c(rep(1, x$input$n_sources))
      # }
      p_means[1:x$input$n_sources,] = diag(x$input$n_sources)
      for(k in 1:x$input$n_obs){
        p_means[k+x$input$n_sources,] = c(colMeans(x$output[[group[i]]]$BUGSoutput$sims.list$p[k,,]))
      }
      colnames(p_means) = x$input$source_names
    }
    size <- c(rep(0.5, x$input$n_sources), rep(0.5, nrow(curr_mix)))
    
    if (ncol(curr_mix) == 1) {
      df <- cbind(data.frame(x = x2, x_lower, x_upper, Source, size, y = Source), p_means)
    } else {
      df <- cbind(data.frame(x = x2, y = y, x_lower, y_lower, x_upper, y_upper, Source, size), p_means)
    }
    
    # Plot for bivariate mixtures
    if (ncol(curr_mix) > 1) {
      if (colour) {
        g <- ggplot(data = df, aes(x = x, y = y, colour = Source)) +
          geom_pie_glyph(slices = colnames(df[,(ncol(df) - x$input$n_sources +1):(ncol(df))]),
                         inherit.aes = TRUE) +
          scale_color_viridis(discrete = TRUE) +
          theme_bw() +
          labs(x = xlab, y = ylab, title = title) +
          geom_errorbarh(aes(xmax = x_upper, xmin = x_lower, height = 0)) +
          # geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
          geom_pointrange(aes(x = x, y = y, ymax = y_upper, ymin = y_lower, shape = Source)) +
          scale_shape_manual(values = 1:nlevels(df$Source)) +
          theme(legend.title = element_blank(), legend.key = element_blank()) +
          guides(color = guide_legend(override.aes = list(linetype = c(rep(0, curr_n_groups), rep(1, x$input$n_sources))))) +
          ggargs  #[(x$input$n_sources+1):(nrow(df)),])
      } else {
        g <- ggplot(data = df, aes(x = x, y = y, colour = Source)) +
          geom_pie_glyph(slices = colnames(df[,(ncol(df) - x$input$n_sources +1):(ncol(df))]),
                         inherit.aes = TRUE) +
          scale_fill_manual(values = 1:nlevels(df$Source)) +
          theme_bw() +
          geom_pie_glyph(slices = colnames(df[,(ncol(df) - x$input$n_sources +1):(ncol(df))])) +
          labs(x = xlab, y = ylab, title = title) +
          geom_errorbarh(aes(xmax = x_upper, xmin = x_lower, height = 0)) +
          # geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
          geom_pointrange(aes(x = x, y = y, ymax = y_upper, ymin = y_lower, shape = Source)) +
          scale_shape_manual(values = 1:nlevels(df$Source)) +
          theme(legend.title = element_blank(), legend.key = element_blank()) +
          guides(color = guide_legend(override.aes = list(linetype = c(rep(0, curr_n_groups), rep(1, x$input$n_sources))))) +
          scale_colour_grey() +
          ggargs
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
          guides(color = guide_legend(override.aes = list(linetype = c(rep(0, curr_n_groups), rep(1, x$n_sources))))) +
          ggargs
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
          guides(color = guide_legend(override.aes = list(linetype = c(rep(0, curr_n_groups), rep(1, x$n_sources))))) +
          ggargs
        
      }
    }
    
    print(g)
    invisible(g)
    }
  }
