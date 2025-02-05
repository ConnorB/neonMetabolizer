#' \code{2stationO2}
#' 
#' Internal function.
#' 
#' @param name description
#' 
#' @returns 
#' 
#' description here
#' 
O2_twoStation <- function(i, oxyup, GPP, z, light, ER, tt, tempup, K600mean, gas,
                          n, osatup, osatdown, lag){
  (oxyup[i] + (
    (GPP/z) *
      ( sum(light[i:(i+lag)]) / sum(light) )
  ) + ER * tt/z +
    (
      Kcor(tempup[i], K600mean, gas = gas, n = n)
    ) * tt * (
      osatup[i] - oxyup[i] + osatdown[i]
    ) / 2
  ) /
    (
      1+ Kcor(tempup[i], K600mean, gas = gas, n = n) * tt / 2
    )
}