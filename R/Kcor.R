#' \code{Kcor} Estimate K, given water temperature and K600 values
#'
#' Internal function. See Wanninkhof 1992 and Hall et al 2016 supplemental material for more detailed information.
#'
#' @param temp Water temperature in degrees C
#' @param K600 Normalized K600 reaeration rate
#'
#' @returns
#'
#' Populate here
#'
#' @references
#'
#' Populate here
#'
#' @example
#'
#' Populate here
#'
#' @export
#'
Kcor <- function(temp, K600) {
  K600 / (600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
