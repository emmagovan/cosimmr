\name{cosimmr-package}
\alias{cosimmr-package}
\alias{cosimmr}
\docType{package}
\title{
  cosimmr: An R package for Stable Isotope Mixing Models
}
\description{
cosimmr is a package that has been developed to allow for running of Stable Isotope Mixing Models in R. It allows for the inclusion of covariates and has been designed to be easy to use for non-expert users. cosimmr uses Fixed Form Variational Bayes to run SIMMs, instead of MCMC. This allows for faster running of models without any issues with convergenge
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Emma Govan, emmagovan@gmail.com.
Andrew C Parnell.

Maintainer: Emma Govan <emmagovan@gmail.com>
}
\references{ add papers
  
}
\keyword{ cosimmr }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
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
#' # Print
#' cosimmr_1
#'
#' # FFVB run
#' cosimmr_1_out <- cosimmr_ffvb(cosimmr_1)
#'
#' # Print it
#' print(cosimmr_1_out)
#'
#' # Summary
#' summary(cosimmr_1_out, type = "correlations")
#' summary(cosimmr_1_out, type = "statistics")
#' ans <- summary(cosimmr_1_out, type = c("quantiles", "statistics"))
#'
#' # Plot
#' plot(cosimmr_1_out, type = "beta_boxplot")
#' plot(cosimmr_1_out, type = "beta_histogram")
#'
  }
}
