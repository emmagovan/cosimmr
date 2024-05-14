#' Alligator Data
#' 
#' Dataset from Nifong et al 2015 which contains 2 food sources, 181 individuals and 2 isotopes. This dataset includes multiple covariates as well as TDF means and sds.
#' 
#' @source <doi:10.1111/1365-2656.12306>
#' 
#' @format A list with the following elements
#' \describe{
#' \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#' \item{ID}{A character vector of unique ID values}
#' \item{tag_ID}{A character vector of tag ID}
#' \item{lat}{A numeric vector of latitude}
#' \item{long}{A numeric vector of longitude}
#' \item{date}{A character vector of date}
#' \item{year}{A numeric vector of year}
#' \item{habitat}{A character vector of habitat}
#' \item{sex}{A character vector of sex}
#' \item{length}{A numeric vector of length in cm}
#' \item{s_class}{A character vector of size class}
#' \item{sex_sclass}{A character vector for sex times size class}
#' \item{source_names}{A character vector of source names}
#' \item{source_means}{A data frame of source means for same tracers as in Mixtures}
#' \item{source_sds}{A data frame of standard deviations of sources for same tracers as in Mixtures}
#' \item{n_sources}{A numeric vector of number of sources}
#' \item{TEF_means}{A data frame of means for TEFs for same tracers as in Mixtures}
#' \item{TEF_sds}{A data frame of sds for TEFs for same tracers as in Mixtures}
#' 
#' }
"alligator_data"

#' Cladocera data from Galloway et al 2015.
#' 
#' Cladocera data from Galloway et al 2015. This dataset has 14 individuals on 7 food sources and 22 tracers. The id column can be used as a covariate. This dataset includes TDFs.
#' 
#' @source <doi:10.1111/fwb.12394>
#' 
#' @format A list with the following elements
#' \describe{
#' \item{id}{numeric vector of ID number}
#' \item{group}{character vector of group ID}
#' \item{mixtures}{Data frame of tracer values. There are 22 fatty acids as tracers in this dataset.}
#' \item{tracer_names}{character vector of tracer names}
#' \item{source_means}{data frame of tracer means for each of the 7 food sources}
#' \item{source_sds}{data frame of tracer sds for each of the 7 food sources}
#' \item{n_sources}{numeric vector of the number of each food source collected}
#' \item{correction_means}{data frame with TDF means for each food source on each tracer}
#' \item{correction_sds}{data frame with TDF sds for each food source on each tracer}
#' 
#' }
"cladocera_data"

#' A smaller version of the Geese stable isotope mixing data set
#'
#' A real Geese data set with 9 observations on 2 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means. Taken from Inger et al (2016). See link for paper
#'
#' @source <doi:10.1111/j.1365-2656.2006.01142.x>
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N, d34S)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'   ...
#'
#' }
"geese_data_day1"

#' Geese stable isotope mixing data set
#'
#' A real Geese data set with 251 observations on 2 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means. Taken from Inger et al (2016). See link for paper
#'
#' @source <doi:10.1111/j.1365-2656.2006.01142.x>
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N, d34S)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'
#' }
"geese_data"

#' Isopod Data
#' 
#' Isopod data from Galloway et al 2014. This dataset has 8 tracers (fatty acids), 30 individuals and 3 food sources. This dataset includes TDFs.
#' 
#' @source <doi:10.3354/meps10860>
#' 
#' @format A list with the following elements
#' \describe{
#' \item{site}{A character vector with name of site for each individual}
#' \item{mixtures}{Data frame with 8 tracer values for 30 individuals}
#' \item{tracer_names}{character vector of tracer names}
#' \item{source_names}{character vector of food source names}
#' \item{source_means}{Data frame of source means with values for each food source on each tracer}
#' \item{source_sds}{Data frame of source sds with values for each food source on each tracer}
#' \item{n_sources}{numeric vector of number of each source obtained}
#' \item{correction_means}{Data frame of TDF means for each food source on each tracer}
#' \item{correction_sds}{Data frame of TDF sds for each food source on each tracer}
#' 
#' }
"iso_data"

#' A simple fake stable isotope mixing data set
#'
#' A simple fake data set with 10 observations on 2 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C and d15N)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'
#' }
"simmr_data_1"

#' A 3-isotope fake stable isotope mixing data set
#'
#' A fake data set with 30 observations on 3 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A three column matrix containing delta 13C, delta 15N, and delta 34S values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N, d34S)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'
#' }
"simmr_data_2"

#' An artificial data set used to indicate effect of priors
#'
#' A fake box data set identified by Fry (2014) as a failing of SIMMs
#' See the link for more interpretation of these data and the output
#'
#' @source <doi:10.3354/meps10535>
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'
#' }
"square_data"
