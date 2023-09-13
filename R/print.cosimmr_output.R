#' Print a simmr output object
#'
#' @param x An object of class \code{cosimmr_output}
#' @param ... Other arguments (not supported)
#'
#' @return Returns a neat summary of the object
#'
#' @seealso  \code{\link{cosimmr_ffvb}} for creating
#' \code{cosimmr_output} objects
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
