#' Calculate two station MCMC during one day. Internal function. Function runs
#' the MCMC and returns posterior distributions
#'
#' O2data    Dataframe containing formatted raw two station data, as
#'           returned by request_NEON in $data
#' upName    Character string denoting name of upstream station (ex. "S1")
#' downName  Character string denoting name of downstream station (ex. "S2")
#' init      Initial guess for the range of GPP and |ER|
#' z         Numeric, average depth on day of modeling
#' tt        Numeric, travel time between stations on day of modeling
#' Kmean     Mean K for informative prior
#' Ksd       Standard deviation of K for informative prior
#' nbatch    Number of MCMC trials
#' scale
twostationpostsum <- function(O2data, upName, downName, start, z, tt, Kmean, Ksd,
                              nbatch, scale) {
  # Create list of unique dates in data
  dateList <- unique(O2data$date)

  # Initialize empty vectors to store results from modeling
  #date <- vector(length = length(unique(data$date)))
  GPP <- vector(mode = "numeric", length = length(unique(O2data$date)))
  GPP.lower <- vector(mode = "numeric", length = length(unique(O2data$date)))
  GPP.upper <- vector(mode = "numeric", length = length(unique(O2data$date)))
  ER <- vector(mode = "numeric", length = length(unique(O2data$date)))
  ER.lower <- vector(mode = "numeric", length = length(unique(O2data$date)))
  ER.upper <- vector(mode = "numeric", length = length(unique(O2data$date)))
  K <- vector(mode = "numeric", length = length(unique(O2data$date)))
  K.lower <- vector(mode = "numeric", length = length(unique(O2data$date)))
  K.upper <- vector(mode = "numeric", length = length(unique(O2data$date)))
  s <- vector(mode = "numeric", length = length(unique(O2data$date)))
  s.lower <- vector(mode = "numeric", length = length(unique(O2data$date)))
  s.upper <- vector(mode = "numeric", length = length(unique(O2data$date)))
  accept <- vector(length = length(unique(O2data$date)))

  # For loop iterating through each day of data
  for (i in 1:length(dateList)){
    # Seperate data into upstream and downstream sections
    updata <- O2data[O2data$horizontalPosition == upName,]
    downdata <- O2data[O2data$horizontalPosition == downName,]

    # NOTE: In Hall et al, the number below was 0.00347222, which corresponds to:
    # 0.00347222 days = 5 minutes, as their O2 sensors took 5-minute readings.
    # Here, our sensors took 15 minute readings. So:
    # 15 minutes = 0.0104166667 days
    # number of 15 min readings bewteen up and down probe corresponding to travel
    # time tt
    lag <- as.numeric(round(tt/0.0104166667))

    # trim the ends of the oxy and temp data by the lag so that oxydown[1]
    # is the value that is the travel time later than oxy up.  The below calls
    # are designed to work with our data structure.

    if (length(updata$DO_mgL) < lag) {
      # If there is less data during a day than the lag interval, move to next day
      message("ERROR: model not computed for ", dateList[i], " as insufficient observations provided.")
      break
    } else{
      # Else continue
      tempup <- updata$WaterTemp_C[1:as.numeric(length(updata$WaterTemp_C)-lag)] # trim the end by the lag
      oxyup <- updata$DO_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
      osat <- updata$DOsat_pct[1:as.numeric(length(updata$WaterTemp_C)-lag)]
      tempdown <- downdata$WaterTemp_C[(1+lag):length(downdata$WaterTemp_C)]
      oxydown <- downdata$DO_mgL[(1+lag):length(downdata$WaterTemp_C)]
      light <- downdata$Light_PAR

      # perform MCMC
      # see documentation on mcmc
      met.post <- mcmc::metrop(tspost, initial = start, nbatch = nbatch,
                               scale = scale,tempup = tempup,
                               tempdown = tempdown, oxyup = oxyup,
                               osat = osat, oxydown = oxydown,  z = z,
                               light = light, tt = tt, Kmean = Kmean,
                               Ksd = Ksd, debug = TRUE)

      # trying to troubleshoot here
      plot(ts(met.post$batch), main = dateList[i])


      # Calculate overall estimates for each day
      gppr <- quantile(met.post$batch[(2000:nbatch),1], c(0.025, 0.5, 0.975))
      err <- quantile(met.post$batch[(2000:nbatch),2], c(0.025, 0.5, 0.975))
      Kr <- quantile(met.post$batch[(2000:nbatch),3], c(0.025, 0.5, 0.975))
      sr <- quantile(met.post$batch[(2000:nbatch),4], c(0.025, 0.5, 0.975))

      # Add results to vectors
      #date[i] <- ymd(unique(data$date)[i])
      GPP[i] <- gppr[2]
      GPP.lower[i] <- gppr[1]
      GPP.upper[i] <- gppr[3]
      ER[i] <- err[2]
      ER.lower[i] <- err[1]
      ER.upper[i] <- err[3]
      K[i] <- Kr[2]
      K.lower[i] <- Kr[1]
      K.upper[i] <- Kr[3]
      s[i] <- sr[2]
      s.lower[i] <- sr[1]
      s.upper[i] <- sr[3]
      accept[i] <- met.post$accept # log likelihood plus priors, should be about 0.2
    }
    i <- i+1
  }
  # Create dataframe of predicted metabolism values
  pred.metab <- data.frame(date = dateList, GPP, GPP.lower, GPP.upper, ER, ER.lower,
                           ER.upper, K, K.lower, K.upper, s, s.lower, s.upper)

  # Call O2TimeSeries function to return modeled O2 values based on median GPP and ER
  # modeling results
  oxymodel <- O2TimeSeries(GPP = pred.metab$GPP, ER = pred.metab$ER,
                           O2data = O2data, Kmean = Kmean, z = z, tt = tt,
                           upName = upName, downName = downName)

  # Create output list to return to user
  output <- list(pred.metab = pred.metab, accept = accept, oxymodel = oxymodel)

  return(output)
}
