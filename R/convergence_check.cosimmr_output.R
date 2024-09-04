#' Creates a plot of mean lower bound to check the FFVB model has finished. Uses
#'  \code{\link{cosimmr_ffvb}} output.
#'
#' Creates diagnostic plot  for an object created with  \code{\link{cosimmr_ffvb}}. 
#'
#' @param object An object of class \code{cosimmr_output} produced by the
#' function \code{\link{cosimmr_ffvb}}
#' @param ...  Not used
#' @return A plot which should show that the mean LB is settling on a value.
#' @author Emma Govan <emmagovan@@gmail.com> Andrew Parnell
#' @seealso See \code{\link{cosimmr_ffvb}}for creating objects suitable for 
#' this function, and many more examples.
#' See also \code{\link{cosimmr_load}} for creating cosimmr objects,
#' \code{\link{plot.cosimmr_input}} for creating isospace plots,
#' \code{\link{plot.cosimmr_output}} for plotting output.
#'
#' @importFrom stats sd cor
#'
#' @examples
#' \donttest{
#' # A simple example with 10 observations, 2 tracers and 4 sources
#'
#' # The data
#' data(geese_data_day1)
#' cosimmr_1 <- with(
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
#'
#' # Plot
#' plot(cosimmr_1)
#'
#'
#' # FFVB run
#' cosimmr_1_out <- cosimmr_ffvb(cosimmr_1)
#'
#' # Summarise
#' convergence_check(cosimmr_1_out)
#' }
#' @export
convergence_check <-
  function(object, ...) {
    lower =(object$output$model$data$ffvb_control$t_W+2)
    upper = object$output$lambdares$iteration
    to_plot = data.frame(Mean_LB = object$output$lambdares$mean_LB_bar[lower:upper], n = 1:length(object$output$lambdares$mean_LB_bar[lower:upper]))

    
    ggplot(data = to_plot, aes(x = n, y = Mean_LB)) + geom_point()

  }


