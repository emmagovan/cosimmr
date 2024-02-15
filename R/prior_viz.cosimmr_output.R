#' Plot the prior distribution for a cosimmr run
#'
#' This function takes the output from \code{\link{cosimmr_ffvb}} and plots the 
#' prior distribution to enable visual
#' inspection. This can be used by itself or together with
#' \code{\link{posterior_predictive}} to visually evaluate the influence of
#' the prior on the posterior distribution.
#'
#' @param cosimmr_out A run of the cosimmr model from \code{\link{cosimmr_ffvb}}
#' @param plot Whether to create a density plot of the prior or not. The simulated prior values are returned as part of the object
#' @param include_posterior Whether to include the posterior distribution on top of the priors. Defaults to TRUE. The posterior returned is of the mean value of covariates
#' @param n_sims The number of simulations from the prior distribution
#' @param scales The type of scale from \code{facet_wrap} allowing for \code{fixed}, \code{free}, \code{free_x}, \code{free_y}
#'
#' @returns A list containing \code{plot}: the ggplot object (useful if requires customisation), and \code{sim}: the simulated prior values which can be compared with the posterior densities
#'
#' @export
#'
#' @examples
#' \donttest{
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
#' # Print
#' simmr_1
#'
#' # MCMC run
#' simmr_1_out <- simmr_mcmc(simmr_1)
#'
#' # Prior predictive
#' prior <- prior_viz(simmr_1_out)
#' head(prior$p_prior_sim)
#' summary(prior$p_prior_sim)
#' }
prior_viz <- function(cosimmr_out,
                      plot = TRUE,
                      include_posterior = TRUE,
                      n_sims = 10000,
                      scales = "free") {
  UseMethod("prior_viz")
}
#' @export
prior_viz.cosimmr_output <- function(cosimmr_out,
                                   plot = TRUE,
                                   include_posterior = TRUE,
                                   n_sims = 10000,
                                   scales = "free") {
    
    # Get group_name
  
    plot_title_1 <- "Prior distributions"
    plot_title_2 <- "Prior and posterior distributions"

    
    # Plot and/or output the prior
    mu_f_mean <- cosimmr_out$output$model$data$mu_f_mean
    sigma_f_sd <- cosimmr_out$output$model$data$sigma_f_sd
    n_sources <- cosimmr_out$input$n_sources
    
    # Now simulate some ps
    p_prior_sim <- matrix(NA, ncol = n_sources, nrow = n_sims)
    for (i in 1:n_sims) {
      f <- rnorm(n_sources, mean = mu_f_mean, sd = sigma_f_sd)
      p_prior_sim[i, ] <- exp(f) / sum(exp(f))
    }
    colnames(p_prior_sim) <- cosimmr_out$input$source_names
    if (plot) {
      
      ############ from plot.output function
      if(cosimmr_out$input$intercept == TRUE){
        x_pred = c(1, rep(0, (ncol(cosimmr_out$input$x_scaled) - 1)))
      } else if(cosimmr_out$input$intercept == FALSE){
        x_pred = c(rep(0, (ncol(x$input$x_scaled))))
      }
      
      thetares= cosimmr_out$output$theta
      K = cosimmr_out$input$n_sources
      n_tracers = cosimmr_out$input$n_tracers
      n_covariates = ncol(cosimmr_out$input$x_scaled)
      n_output = ncol(cosimmr_out$output$theta)
      
      
      
      sigma <- cosimmr_out$output$BUGSoutput$sims.list$sigma#(sqrt(exp(thetares[,(K*n_covariates + 1):(K*n_covariates + n_tracers)]))
      
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
      
      colnames(p_sample) = cosimmr_out$input$source_names
      
      df_p_mean <- reshape2::melt(p_sample)
      
      
      ##############################
      
      df <- reshape2::melt(p_prior_sim)
      colnames(df) <- c("Num", "Source", "Proportion")
      df$Type <- "Prior"
     # out_all <- cosimmr_out$output$BUGSoutput$sims.list$p
      df2 <- df_p_mean
      colnames(df2) <- c("Num", "Source", "Proportion")
      df2$Type <- "Posterior"
      df_all <- rbind(df2, df)
      if (include_posterior) {
        g <- ggplot(
          df_all,
          aes(
            x = Proportion,
            y = after_stat(density),
            fill = Source,
            linetype = Type
          )
        ) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(alpha = 0.5) +
          theme_bw() +
          ggtitle(plot_title_2) +
          ylab("Density") +
          facet_wrap("~ Source", scales = scales)
      } else {
        g <- ggplot(
          df,
          aes(
            x = Proportion,
            y = after_stat(density),
            fill = Source
          )
        ) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(alpha = 0.5, linetype = 0) +
          theme_bw() +
          ggtitle(plot_title_1) +
          ylab("Density") +
          facet_wrap("~ Source", scales = scales)
      }
      print(g)
    }
    
    # Return the simulations
    if (exists("g")) {
      invisible(list(plot = g, p_prior_sim = p_prior_sim))
    } else {
      invisible(p_prior_sim)
    }
  
}
