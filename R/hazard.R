##############################################################
#' Constant hazard function
#'
#' Hazard function for use with gaussian underlying
#' distribution.
#'
#' @param lambda The parameter for the hazard function.
#'
#' @param r The current R vector length.
#'
#' @return A vector of the hazard function for the length
#' of the current R vector.
#'
#' @docType methods
#'
#' @export
#'
#' @examples
#' H<- const_hazard(10, 1/100)
##############################################################
const_hazard <-
function(r, lambda) {
  return(1/lambda * rep(1,r))
}

