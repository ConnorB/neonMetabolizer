#' \code{O2TimeSeries} Return modeled oxygen time series from median GPP, ER estimates
#'
#' Internal function.
#'
#' @param GPP Description
#' @param ER Description
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
O2TimeSeries <- function(GPP, ER, timeup, timedown, oxyup, z, light,
                         tt, tempup, K600mean, gas, n, osatup, osatdown,
                         lag, oxydown) {
  # Initialize an empty vector
  modeledO2 <- numeric(length(oxyup))

  # Calculate metabolism at each timestep
  for (i in 1:length(oxyup)) {
    # Check if non-NA data on date, if GPP is NA, go to next date in sequence
    if(is.na(GPP)){
      # If there is no GPP estimate, move on to next row in dataframe
      next
    } else{
      modeledO2[i] <- O2_twoStation(i = i, oxyup = oxyup, GPP = GPP, z = z,
                                    light = light, ER = ER, tt = tt,
                                    tempup = tempup, K600mean = K600mean,
                                    gas = gas, n = n, osatup = osatup,
                                    osatdown = osatdown, lag = lag)
    } # close else
  }

  # Convert time from chron to posixct
  timeup <- as.POSIXlt(timeup, origin = "1970-01-01")
  timedown <- as.POSIXlt(timedown, origin = "1970-01-01")

  gasmodel <- data.frame(timeup, timedown,
                         oxydown, oxyup, modeledO2)
  return(gasmodel)
}
