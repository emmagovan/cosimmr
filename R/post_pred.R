#' Plot the posterior predictive distribution for a cosimmr run
#'
#' This function takes the output from  \code{\link{cosimmr_ffvb}} and plots 
#' the posterior predictive distribution to enable visualisation of model fit.
#' The simulated posterior predicted values are returned as part of the object 
#' and can be saved for external use
#'
#' @param cosimmr_out A run of the cosimmr model from \code{\link{cosimmr_ffvb}}.
#' @param prob The probability interval for the posterior predictives. The default is 0.5 (i.e. 50pc intervals)
#' @param plot_ppc Whether to create a bayesplot of the posterior predictive or not.
#' @param ind The individual number you wish to plot. Defaults to 1.
#' @param n_samples The number of samples you wish to generate for y_pred. Defaults to 3600.
#'
#'@return plot of posterior predictives and simulated values
#'
#'#' @author Emma Govan <emmagovan@@gmail.com> Andrew Parnell
#'
#' @seealso \code{\link{cosimmr_ffvb}} for creating objects suitable for this
#' function
#' 
#' @importFrom bayesplot ppc_intervals
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(geese_data_day1)
#' cosimmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ c(1,2,3,2,1,2,3,2,1),
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
#' # Print
#' cosimmr_1
#'
#' # FFVB run
#' cosimmr_1_out <- cosimmr_ffvb(cosimmr_1)
#'
#' # Prior predictive
#' post_pred <- posterior_predictive(cosimmr_1_out)
#' }
posterior_predictive <- function(cosimmr_out,
                                 prob = 0.5, 
                                 ind = 1,
                                 plot_ppc = TRUE,
                                 n_samples = 3600) {
  UseMethod("posterior_predictive")
}
#' @export
posterior_predictive.cosimmr_output <- function(cosimmr_out,
                                              prob = 0.5,
                                              ind = 1,
                                              plot_ppc = TRUE,
                                              n_samples = 3600) {
 
  #Okay so what I want to do here is to  take theta matrix and simulate from likelihood
  #So we have y ~ N(mu_y, sigma_y) and we want to sample from that?
  #And then we take our actual original y data and compare? I think
  K = cosimmr_out$input$n_sources
  n_covariates = ncol(cosimmr_out$input$x_scaled)
  n_tracers = cosimmr_out$input$n_tracers
  theta = cosimmr_out$output$theta
  n_output = nrow(theta)
  x_pred = cosimmr_out$input$x_scaled
  n_obs = cosimmr_out$input$n_obs
  
#Need to specify an individual I guess?
  p = cosimmr_out$output$BUGSoutput$sims.list$p[ind,,]
  sigma = (cosimmr_out$output$BUGSoutput$sims.list$sigma)
  #Need to add proper loops here to do it properly
  #But think the general idea is then do rnorm(mean_y, sigma_y)
  
  q = cosimmr_out$input$concentration_means
  mu_s = cosimmr_out$input$source_means
  sigma_s = cosimmr_out$input$source_sds
  mu_c = cosimmr_out$input$correction_means
  sigma_c = cosimmr_out$input$correction_sds
  
  mean = matrix(NA, nrow = n_obs, ncol = n_tracers)
  sd = matrix(NA, nrow = n_obs, ncol = n_tracers)
  for (i in 1:n_obs) {
    for (j in 1:n_tracers) {
    mean[i,j] = sum(p[i, ] * q[, j] * (mu_s[, j] + mu_c[, j])) /
          sum(p[i, ] * q[, j])
    sd[i,j] = sqrt(sum(p[i, ]^2 * q[, j]^2 * (sigma_s[, j]^2 + sigma_c[, j]^2)) /
          sum(p[i, ]^2 * q[, j]^2) + sigma[i,j])
      

      
    }
  }

  
  
  y_post_pred = array(NA, dim = c(n_samples, n_obs, n_tracers))
  for(k in 1:n_tracers){
    for(i in 1:n_obs){
y_post_pred[,i,k] = stats::rnorm(n_samples, mean = mean[i,k], sd = sd[i,k])
    }
  }
  


  
  # Make is look nicer
  low_prob <- 0.5 - prob / 2
  high_prob <- 0.5 + prob / 2
   y_post_pred_ci <- apply(y_post_pred,
                           2:3,
                           "quantile",
                           prob = c(low_prob, high_prob)
   )
  y_post_pred_out <- data.frame(
    interval = matrix(y_post_pred_ci,
                      ncol = cosimmr_out$input$n_tracers,
                      byrow = TRUE
    ),
    data = as.vector(cosimmr_out$input$mixtures)
  )
  
  y_post_pred_out$outside <- y_post_pred_out[, 3] > y_post_pred_out[, 2] |
    y_post_pred_out[, 3] < y_post_pred_out[, 1]
  prop_outside <- mean(y_post_pred_out$outside)
  
  if (plot_ppc) {
    y_rep <- y_post_pred
    dim(y_rep) <- c(dim(y_post_pred)[1], dim(y_post_pred)[2] * dim(y_post_pred)[3])
 
    curr_mix <- cosimmr_out$input$mixtures
    bayesplot::color_scheme_set("viridis")
    bayesplot::bayesplot_theme_set(new = theme_bw())
    g <- ppc_intervals(
      y = unlist(as.vector(curr_mix)),
      yrep = y_rep,
      x = rep(1:nrow(curr_mix), cosimmr_out$input$n_tracers),
      prob = prob,
      fatten = 1
    ) + ggplot2::ylab("Tracer value") +
      ggplot2::xlab("Observation") +
      ggplot2::ggtitle(paste0(prob * 100, "% posterior predictive")) +
      ggplot2::scale_x_continuous(breaks = 1:cosimmr_out$input$n_obs)
    print(g)
  }
  # Return the simulations
  invisible(list(
    table = y_post_pred_out,
    prop_outside = prop_outside
  ))
}
