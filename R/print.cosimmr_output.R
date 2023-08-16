#' Print a simmr output object
#'
#' @param x An object of class \code{simmr_output}
#' @param ... Other arguments (not supported)
#'
#' @return Returns a neat summary of the object
#'
#' @seealso \code{\link{simmr_mcmc}} and \code{\link{simmr_ffvb}} for creating
#' \code{simmr_output} objects
#' @export
print.cosimmr_output <-
  function(x, ...) {
    if (inherits(x, "cosimmr_output") == TRUE) {
         if (inherits(x, "ffvb") == TRUE) {
        print(x$input)
        message("The input data has been run via simmr_ffvb and has produced ")
        message(nrow(x$output$BUGSoutput$sims.list$sigma), " samples.")
      }
    }
  }
