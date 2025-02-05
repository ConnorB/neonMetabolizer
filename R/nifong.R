#' \code{NifongEqn} Model N2 using two-station equation from Nifong et al. 2020
#'
#' Internal function.
#'
#' @param name description
#'
#' @returns
#'
#' description here
#'
nifong <- function(i, n2up, NConsume, z, light, DN, tt, tempup, K600mean, gas,
                   n, nsatup, nsatdown, lag){
  (n2up[i] + (
    (NConsume/z) *
      ( sum(light[i:(i+lag)]) / sum(light) )
  ) + DN * tt/z +
    (
      Kcor(tempup[i], K600mean, gas = gas, n = n)
    ) * tt * (
      nsatup[i] - n2up[i] + nsatdown[i]
    ) / 2
  ) /
    (
      1 + Kcor(tempup[i], K600mean, gas = gas, n = n) * tt / 2
    )
}
