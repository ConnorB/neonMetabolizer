#' \code{lightIndependent}
#'
#' Internal function.
#'
#' @param name description
#'
#' @returns
#'
#' description here
#'
lightIndependent <- function(i, n2up, NConsume, z, DN, tt, tempup,
                             K600mean, gas, n, nsatup, nsatdown){
  (n2up[i] +
     NConsume * tt/z +
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
