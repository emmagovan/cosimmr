#' Plot the posterior predictive distribution for a cosimmr run
#'
#' This function takes the output from  \code{\link{cosimmr_ffvb}} and plots 
#' the posterior predictive distribution to enable visualisation of model fit.
#' The simulated posterior predicted values are returned as part of the object 
#' and can be saved for external use
#'
#' @param cosimmr_out A run of the simmr model from \code{\link{cosimmr_ffvb}}.
#' @param covariate Which covariate to run it for
#' @param prob The probability interval for the posterior predictives. The default is 0.5 (i.e. 50pc intervals)
#' @param plot_ppc Whether to create a bayesplot of the posterior predictive or not.
#'
#'@return plot of posterior predictives and simulated values

#' @importFrom bayesplot ppc_intervals
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(geese_data_day1)
#' simmr_1 <- with(
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
#' plot(simmr_1)
#'
#' # Print
#' simmr_1
#'
#' # FFVB run
#' simmr_1_out <- cosimmr_ffvb(simmr_1)
#'
#' # Prior predictive
#' post_pred <- posterior_predictive(simmr_1_out)
#' }
posterior_predictive <- function(simmr_out,
                                 group = 1,
                                 prob = 0.5,
                                 plot_ppc = TRUE) {
  UseMethod("posterior_predictive")
}
#' @export
posterior_predictive.simmr_output <- function(simmr_out,
                                              group = 1,
                                              prob = 0.5,
                                              plot_ppc = TRUE) {
 
  #Okay so what I want to do here is to  take theta matrix and simulate from likelihood
  #So we have y ~ N(mu_y, sigma_y) and we want to sample from that?
  #And then we take our actual original y data and compare? I think
  K = simmr_out$input$n_sources
  n_covariates = ncol(simmr_out$input$x_scaled)
  n_tracers = simmr_out$input$n_tracers
  theta = simmr_out$output$theta
  n_output = nrow(theta)
  x_pred = simmr_out$input$x_scaled
  n_obs = simmr_out$input$n_obs
  
  sigma <- (1/sqrt(theta[,(K*n_covariates + 1):(K*n_covariates + n_tracers)]))
  beta = array(theta[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))
  
  #So I think what we want to do here is to generate p using beta? And then
  #Put all the conc means etc together to create y
  #And then sample from y
  
  f <- array(NA, dim = c(nrow(x_pred), K, n_output))
  p = array(NA, dim =  c(n_obs, n_output, K))
  
  for(s in 1:n_output){
    f[,,s] = as.matrix(x_pred) %*% beta[s,,]
  }
  
  for(j in 1:n_output){
    for (k in 1:n_obs) {
      p[k,j, ] <- exp(f[k,1:K, j]) / (sum((exp(f[k,1:K, j]))))
    }
  }
  
  #Need to add proper loops here to do it properly
  #But think the general idea is then do rnorm(mean_y, sigma_y)
  
  mean_y = sum((p * q) * (mu_s + mu_c)) / sum(p*q)
  sigma_y = sum((p^2 * q^2) * (sigma_s^2 + sigma_c)^2) / sum(p^2*q^2)
  

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
                      ncol = simmr_out$input$n_tracers,
                      byrow = TRUE
    ),
    data = as.vector(simmr_out$input$mixtures[simmr_out$input$group_int == group, ])
  )
  
  y_post_pred_out$outside <- y_post_pred_out[, 3] > y_post_pred_out[, 2] |
    y_post_pred_out[, 3] < y_post_pred_out[, 1]
  prop_outside <- mean(y_post_pred_out$outside)
  
  if (plot_ppc) {
    y_rep <- y_post_pred
    dim(y_rep) <- c(dim(y_post_pred)[1], dim(y_post_pred)[2] * dim(y_post_pred)[3])
    curr_rows <- which(simmr_out$input$group_int == group)
    curr_mix <- simmr_out$input$mixtures[curr_rows, , drop = FALSE]
    bayesplot::color_scheme_set("viridis")
    bayesplot::bayesplot_theme_set(new = theme_bw())
    g <- ppc_intervals(
      y = unlist(as.vector(curr_mix)),
      yrep = y_rep,
      x = rep(1:nrow(curr_mix), simmr_out$input$n_tracers),
      prob = prob,
      fatten = 1
    ) + ggplot2::ylab("Tracer value") +
      ggplot2::xlab("Observation") +
      ggplot2::ggtitle(paste0(prob * 100, "% posterior predictive")) +
      ggplot2::scale_x_continuous(breaks = 1:simmr_out$input$n_obs)
    print(g)
  }
  # Return the simulations
  invisible(list(
    table = y_post_pred_out,
    prop_outside = prop_outside
  ))
}
