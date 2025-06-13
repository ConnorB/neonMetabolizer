#' \code{BlendedEqn2}
#'
#' Internal function.
#'
#' @param name description
#'
#' @returns
#'
#' description here
#'
blended2 <- function(i, n2up, NConsume, NFix, z, light, DN, tt, tempup, K600mean, gas,
                     n, nsatup, nsatdown, lag){
  (n2up[i] +
     (NConsume - NFix) * tt/z +
     (
       (NFix / z) *
         (sum(light[i:(i+lag)]) / sum(light) )
     ) +
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
