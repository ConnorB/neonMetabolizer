#' \code{O2TimeSeries} Return modeled oxygen time series from median GPP, ER estimates
#'
#' Internal function.
#'
#' @param GPP Description
#' @param ER Description
#' @param O2data Dataframe of cleaned raw two station data (ex. "TS_S1S2")
#' @param Kmean Description
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
O2TimeSeries <- function(GPP, ER, O2data, Kmean, z, tt, upName, downName) {
  # Ungroup O2data
  O2data <- O2data %>% dplyr::ungroup()

  #number of 15 min readings bewteen up and down probe corresponding
  # to travel time tt
  lag <- as.numeric(round(tt/0.0104166667))

  # trim the ends of the oxy and temp data by the lag so that oxydown[1]
  # is the value that is the travel time later than oxy up. The below
  # calls are designed to work with our data structure.

  # Seperate data into upstream and downstream sections
  updata <- O2data[O2data$horizontalPosition == upName,]
  downdata <- O2data[O2data$horizontalPosition == downName,]

  tempup <- updata$WaterTemp_C[1:as.numeric(length(updata$WaterTemp_C)-lag)] # trim the end by the lag
  tempdown <- downdata$WaterTemp_C[(1+lag):length(downdata$WaterTemp_C)]

  oxyup <- updata$DO_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
  # define osat
  osat <- updata$DOsat_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
  oxydown <- downdata$DO_mgL[(1+lag):length(downdata$WaterTemp_C)]

  timeup <- updata$dtime[1:(length(updata$WaterTemp_C)-lag)]
  timedown <- downdata$dtime[(1+lag):length(downdata$WaterTemp_C)]

  light <- downdata$Light_PAR

  # Initialize an empty vector
  modeledO2 <- numeric(length(oxyup))
  # Calculate metabolism at each timestep
  for (i in 1:length(oxyup)) {
    modeledO2[i] <- (oxyup[i] + ((GPP/z)*(sum(light[i:(i+lag)]) / sum(light))) +
                       ER*tt/z +
                       (Kcor(tempup[i],Kmean))*tt*(osat[i] - oxyup[i] +  osat[i])/2) /
      (1 + Kcor(tempup[i],Kmean)*tt/2)
  }

  # Convert time from chron to posixct
  timeup <- as.POSIXlt(timeup, origin = "1970-01-01")
  timedown <- as.POSIXlt(timedown, origin = "1970-01-01")

  oxymodel <- data.frame(timeup, timedown, oxydown, oxyup, modeledO2)
  return(oxymodel)
}
