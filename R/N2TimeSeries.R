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
N2TimeSeries <- function(NConsume, DN, data, K600mean, z, tt, upName, downName, gas, n) {
  # Ungroup data
  data <- data %>% dplyr::ungroup()

  #number of 15 min readings bewteen up and down probe corresponding
  # to travel time tt
  lag <- as.numeric(round(tt/0.0104166667))

  # trim the ends of the oxy and temp data by the lag so that oxydown[1]
  # is the value that is the travel time later than oxy up. The below
  # calls are designed to work with our data structure.

  # Seperate data into upstream and downstream sections
  updata <- data[data$horizontalPosition == upName,]
  downdata <- data[data$horizontalPosition == downName,]

  tempup <- updata$WaterTemp_C[1:as.numeric(length(updata$WaterTemp_C)-lag)] # trim the end by the lag
  tempdown <- downdata$WaterTemp_C[(1+lag):length(downdata$WaterTemp_C)]

  n2up <- updata$N2_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
  # define osat
  nsat <- updata$N2sat_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
  n2down <- downdata$N2_mgL[(1+lag):length(downdata$WaterTemp_C)]

  timeup <- updata$solarTime[1:(length(updata$WaterTemp_C)-lag)]
  timedown <- downdata$solarTime[(1+lag):length(downdata$WaterTemp_C)]

  light <- downdata$Light_PAR

  # Initialize an empty vector
  modeledN2 <- numeric(length(n2up))
  # Calculate metabolism at each timestep
  for (i in 1:length(n2up)) {
    # Check if non-NA data on date, if GPP is NA, go to next date in sequence
    if(is.na(NConsume)){
      # If there is no GPP estimate, move on to next row in dataframe
      next
    } else{
      modeledN2[i] <- (n2up[i] + ((NConsume/z)*(sum(light[i:(i+lag)]) / sum(light))) +
                         DN*tt/z +
                         (Kcor(tempup[i],K600mean, gas = gas, n = n))*tt*(nsat[i] - n2up[i] +  nsat[i])/2) /
        (1 + Kcor(tempup[i],K600mean, gas = gas, n = n)*tt/2)
    } # close else
  }

  # Convert time from chron to posixct
  timeup <- as.POSIXlt(timeup, origin = "1970-01-01")
  timedown <- as.POSIXlt(timedown, origin = "1970-01-01")

  gasmodel <- data.frame(timeup, timedown,
                         n2down, n2up, modeledN2)
  return(gasmodel)
}
