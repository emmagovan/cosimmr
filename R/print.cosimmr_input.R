#' Print simmr input object
#'
#' @param x An object of class \code{simmr_input}
#' @param ... Other arguments (not supported)
#'
#' @return A neat presentation of your simmr object.
#' @export print.cosimmr_input
print.cosimmr_input <-
  function(x, ...) {
    message("This is a valid simmr input object with ")
    message(paste(x$n_obs, " observations, "))
    message(paste(ncol(x$x_scaled), " covariates, "))
    message(paste(x$n_tracers, "tracers, and "))
    message(paste(x$n_sources, "sources.\n"))
    message(" The formula is ") 
    message( "c(",paste(colnames(x$mixtures), collapse = ", "), ")~",
            paste(colnames(x$x_scaled), collapse = " + "))
    message("\nThe source names are: ")
    print(x$source_names, sep = ", ")
    message(".\n")
    message("The tracer names are: ")
    print(colnames(x$mixtures), sep = ", ")
    message("\n\n")
  }
