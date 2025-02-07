#' \code{ReisingerEqn} Model N2 using two-station equation from reisinger et al 2016
#'
#' Internal function.
#'
#' @param name description
#'
#' @returns
#'
#' description here
#'
reisinger <- function(i, n2up, z, DN, tt, tempup, K600mean, gas,
                   n, nsatup, nsatdown, lag){
  (n2up[i] +
     DN * tt/z +
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
