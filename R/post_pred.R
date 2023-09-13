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
posterior_predictive <- function(cosimmr_out,
                                 prob = 0.5,
                                 plot_ppc = TRUE,
                                 n_samples = 3600) {
  UseMethod("posterior_predictive")
}
#' @export
posterior_predictive.cosimmr_output <- function(cosimmr_out,
                                              group = 1,
                                              prob = 0.5,
                                              plot_ppc = TRUE) {
 
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
  
  tau <- theta[,(K*n_covariates + 1):(K*n_covariates + n_tracers)]
  beta_all <- theta[,1:(n_covariates * K)]
  
 # beta_all = array(theta[,1:(n_covariates * K)], dim = c(n_output, n_covariates, K))
  
  beta_average = colMeans(beta_all)
  beta_average_mat = matrix(beta_average, nrow = n_covariates, byrow = TRUE)
  tau_average = colMeans(tau)
  #So I think what we want to do here is to generate p using beta? And then
  #Put all the conc means etc together to create y
  #And then sample from y
  
  #f <- array(NA, dim = c(nrow(x_pred), K, n_output))
  #p = array(NA, dim =  c(n_obs, n_output, K))
  f = matrix(NA, nrow = nrow(x_pred), ncol = K)
  p = matrix(NA, nrow = n_obs, ncol = K)
  
  #Probably want to get the average of p over n_output?
  
  f = (x_pred) %*% beta_average_mat
  
  
  
  
  # for(s in 1:n_output){
  #  f[,,s] = as.matrix(x_pred) %*% beta[s,,]
  #   
  # }
  
  for(j in 1:n_obs){
   # for (k in 1:n_obs) {
      p[j, ] <- exp(f[j,]) / (sum((exp(f[j,]))))
  #  }
  }
  
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
          sum(p[i, ]^2 * q[, j]^2) + 1/tau[j])
      

      
    }
  }

  
  
  y_post_pred = array(NA, dim = c(n_samples, n_obs, n_tracers))
  for(k in 1:n_tracers){
    for(i in 1:n_obs){
y_post_pred[,i,k] = rnorm(n_samples, mean = mean[i,k], sd = sd[i,k])
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
