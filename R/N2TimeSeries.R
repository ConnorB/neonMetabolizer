#' \code{N2TimeSeries} Return modeled oxygen time series from median GPP, ER estimates
#'
#' Internal function.
#'
#' @param NConsume Description
#' @param DN Description
#' @param data Dataframe of cleaned raw two station data (ex. "TS_S1S2")
#' @param K600mean Description
#' @param z Description
#' @param tt Description
#' @param upName Description
#' @param downName Description
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
N2TimeSeries <- function(NConsume = NULL, DN, timeup, timedown,
                         NFix = NULL, NOther = NULL,
                         n2up, z, light, tt, tempup, K600mean, gas, n,
                         nsatup, nsatdown, lag, n2down, eqn) {

  # Initialize an empty vector
  modeledN2 <- numeric(length(n2up))

  # Calculate N2 concentrations at each timestep
  for (i in 1:length(n2up)) {
    # Check if non-NA data on date, if GPP is NA, go to next date in sequence
    if(is.na(DN)){
      # If there is no GPP estimate, move on to next row in dataframe
      next
    } else{
      if(eqn == "Nifong_et_al_2020"){
        for (i in 1:length(n2up)){
          # (this function from nifong et al 2020)
          modeledN2[i] <- nifong(i = i, n2up = n2up, NConsume = NConsume, z = z,
                                 light = light, DN = DN, tt = tt, tempup = tempup,
                                 K600mean = K600mean, gas = gas, n = n,
                                 nsatup = nsatup, nsatdown = nsatdown, lag = lag)
        }
      }
      if(eqn == "light_independent"){
        for (i in 1:length(n2up)){
          modeledN2[i] <- lightIndependent(i = i, n2up = n2up,
                                           NConsume = NConsume, z = z,
                                           DN = DN, tt = tt, tempup = tempup,
                                           K600mean = K600mean, gas = gas,
                                           n = n, nsatup = nsatup,
                                           nsatdown = nsatdown)
        }
      }
      if(eqn == "Reisinger_et_al_2016"){
        for (i in 1:length(n2up)){
          modeledN2[i] <- reisinger(i = i, n2up = n2up, z = z,
                                           DN = DN, tt = tt, tempup = tempup,
                                           K600mean = K600mean, gas = gas,
                                           n = n, nsatup = nsatup,
                                           nsatdown = nsatdown)
        }
      }
      if(eqn %in% c("blended1", "blended3")){
        for (i in 1:length(n2up)){
          modeledN2[i] <- blended1(i = i, n2up = n2up, light = light,
                                           NOther = NOther, NFix = NFix, z = z,
                                           DN = DN, tt = tt, tempup = tempup,
                                           K600mean = K600mean, gas = gas,
                                           n = n, nsatup = nsatup,
                                           nsatdown = nsatdown, lag = lag)
        }
      }
      if(eqn == "blended2"){
        for (i in 1:length(n2up)){
          modeledN2[i] <- blended2(i = i, n2up = n2up, light = light,
                                   NOther = NOther, NFix = NFix, z = z,
                                   DN = DN, tt = tt, tempup = tempup,
                                   K600mean = K600mean, gas = gas,
                                   n = n, nsatup = nsatup,
                                   nsatdown = nsatdown, lag = lag)
        }
      }
    } # close else
  }

  # Convert time from chron to posixct
  timeup <- as.POSIXlt(timeup, origin = "1970-01-01")
  timedown <- as.POSIXlt(timedown, origin = "1970-01-01")

  gasmodel <- data.frame(timeup, timedown,
                         n2down, n2up, modeledN2)
  return(gasmodel)
}
