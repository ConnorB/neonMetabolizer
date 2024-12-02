#' \code{Kcor} Estimate kO2, given water temperature and K600 values
#'
#' Internal function. See Wanninkhof 1992 and Hall et al 2016 supplemental material for more detailed information.
#'
#' @param temp Water temperature in degrees C
#' @param K600 Normalized K600 reaeration rate
#' @param gas Either "O2" or "N2", used for calculating Schmidt number
#' @param n Default is 0.5 for wavy water, but 0.667 refers to still water, see Jahne 1987
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
Kcor <- function(temp, K600, gas, n = 0.5) {
  # Values are from metaanalysis in Raymond et al 2012, Table 1
  # Original Hall et al 2016 script used values from Wilke and Chang 1955
  if(gas == "N2"){
    A <- 1615
    B <- -92.15
    C <- 2.349
    D <- -0.0240
  }
  if(gas == "O2"){
    A <- 1568
    B <- -86.04
    C <- 2.142
    D <- -0.0216
  }
  # Calculate Schmidt number for gas
  sc_gas <- A + B*temp + C*(temp^2) + D*(temp^3)

  # Convert K600 values into k_gas values
  k_gas <- K600 / ((600/sc_gas)^(-n))
  return(k_gas)
}
