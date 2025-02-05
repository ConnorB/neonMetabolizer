#' \code{twostationpostsum} Calculate two station MCMC or Bayes during one day.
#'
#' Internal function. Function runs the MCMC and returns posterior distributions
#'
#' @param data Dataframe containing formatted raw two station data, as returned by request_NEON in $data
#' @param upName Character string denoting name of upstream station (ex. "S1")
#' @param downName Character string denoting name of downstream station (ex. "S2")
#' @param start Description
#' @param z Numeric, average depth on day of modeling
#' @param tt Numeric, travel time between stations on day of modeling
#' @param K600mean Mean K600 for informative prior
#' @param K600sd Standard deviation of K600 for informative prior
#' @param nbatch Number of MCMC trials
#' @param scale Description
#' @param modType Description
#' @param gas Either "O2" or "N2" denoting gas type
#'
#' @returns Populate here
#'
#' @references Populate here
#'
#' @example Populate here
twostationpostsum <- function(data, upName, downName, start, z, tt, K600mean, K600sd,
                              nbatch, scale, modType, gas, n, eqn) {
  # Create list of unique dates in data
  dateList <- unique(data$date)

  # Initialize empty vectors to store results from modeling
  #date <- vector(length = length(unique(data$date)))
  K600 <- vector(mode = "numeric", length = length(unique(data$date)))
  K600.lower <- vector(mode = "numeric", length = length(unique(data$date)))
  K600.upper <- vector(mode = "numeric", length = length(unique(data$date)))
  s <- vector(mode = "numeric", length = length(unique(data$date)))
  s.lower <- vector(mode = "numeric", length = length(unique(data$date)))
  s.upper <- vector(mode = "numeric", length = length(unique(data$date)))
  accept <- vector(length = length(unique(data$date)))

  if(gas == "O2"){
    GPP <- vector(mode = "numeric", length = length(unique(data$date)))
    GPP.lower <- vector(mode = "numeric", length = length(unique(data$date)))
    GPP.upper <- vector(mode = "numeric", length = length(unique(data$date)))
    ER <- vector(mode = "numeric", length = length(unique(data$date)))
    ER.lower <- vector(mode = "numeric", length = length(unique(data$date)))
    ER.upper <- vector(mode = "numeric", length = length(unique(data$date)))
  }
  if(gas == "N2"){
    NConsume <- vector(mode = "numeric", length = length(unique(data$date)))
    NConsume.lower <- vector(mode = "numeric", length = length(unique(data$date)))
    NConsume.upper <- vector(mode = "numeric", length = length(unique(data$date)))
    DN <- vector(mode = "numeric", length = length(unique(data$date)))
    DN.lower <- vector(mode = "numeric", length = length(unique(data$date)))
    DN.upper <- vector(mode = "numeric", length = length(unique(data$date)))

    if(grepl(pattern = "blende.+", eqn)){
      NOther <- vector(mode = "numeric", length = length(unique(data$date)))
      NOther.lower <- vector(mode = "numeric", length = length(unique(data$date)))
      NOther.upper <- vector(mode = "numeric", length = length(unique(data$date)))
      NFix <- vector(mode = "numeric", length = length(unique(data$date)))
      NFix.lower <- vector(mode = "numeric", length = length(unique(data$date)))
      NFix.upper <- vector(mode = "numeric", length = length(unique(data$date)))
    }
  }

  i <- 1
  # For loop iterating through each day of data
  for (i in 1:length(dateList)){
    # Seperate data into upstream and downstream sections
    updata <- data[data$horizontalPosition == upName,]
    downdata <- data[data$horizontalPosition == downName,]

    # Get elapsed time between sensor readings in units of days
    tdiff <- as.numeric(difftime(updata$dateTime_local[2],
                                 updata$dateTime_local[1],
                                 units = "days"))

    # Divide travel time by sensor timestep, then round to get timestep lag
    lag <- as.numeric(round(tt/tdiff))

    # In the event that travel time is impossibly fast between stations
    # (lag = 0) set lag = 1
    if(lag == 0){
      lag <- 1
    }

    # trim the ends of the oxy and temp data by the lag so that oxydown[1]
    # is the value that is the travel time later than oxy up.  The below calls
    # are designed to work with our data structure.
    if(length(updata$DO_mgL) < lag) {
      # If there is less data during a day than the lag interval, move to next day
      message("ERROR: model not computed for ",
              dateList[i], " as insufficient observations provided.")
      # Add NA values to dataframe
      K600[i] <- NA
      K600.lower[i] <- NA
      K600.upper[i] <- NA
      s[i] <- NA
      s.lower[i] <- NA
      s.upper[i] <- NA
      accept[i] <- NA
      if(gas=="O2"){
        GPP[i] <- NA
        GPP.lower[i] <- NA
        GPP.upper[i] <- NA
        ER[i] <- NA
        ER.lower[i] <- NA
        ER.upper[i] <- NA
      }
      if(gas=="N2"){
        NConsume[i] <- NA
        NConsume.lower[i] <- NA
        NConsume.upper[i] <- NA
        DN[i] <- NA
        DN.lower[i] <- NA
        DN.upper[i] <- NA
        if(grepl(pattern = "blende.+", eqn)){
          NOther[i] <- NA
          NOther.lower[i] <- NA
          NOther.upper[i] <- NA
          NFix[i] <- NA
          NFix.lower[i] <- NA
          NFix.upper[i] <- NA
        }
      }
      next
      } else{
        # Else continue
        timeup <- updata$solarTime[1:as.numeric(length(updata$WaterTemp_C)-lag)]
        timedown <- downdata$solarTime[(1+lag):length(downdata$WaterTemp_C)]
        tempup <- updata$WaterTemp_C[1:as.numeric(length(updata$WaterTemp_C)-lag)] # trim the end by the lag
        tempdown <- downdata$WaterTemp_C[(1+lag):length(downdata$WaterTemp_C)]
        lightup <- updata$Light_PAR[1:as.numeric(length(updata$Light_PAR)-lag)]
        lightdown <- downdata$Light_PAR[(1+lag):length(downdata$Light_PAR)]
        lightmean <- (lightup + lightdown)/2
        # Script needs to sum i:i+lag, problem for end of time string.
        finallight <- lightmean[length(lightmean)]
        # Extrapolate final light as last light values
        lightmean <- c(lightmean, rep(finallight, 1+lag))
        light <- lightmean

        #light <- lightdown
        if(gas == "O2"){
          oxyup <- updata$DO_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          osatup <- updata$DOsat_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          oxydown <- downdata$DO_mgL[(1+lag):length(downdata$WaterTemp_C)]
          osatdown <- downdata$DOsat_mgL[(1+lag):length(downdata$WaterTemp_C)]
        }
        if(gas == "N2"){
          n2up <- updata$N2_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          nsatup <- updata$N2sat_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          n2down <- downdata$N2_mgL[(1+lag):length(downdata$WaterTemp_C)]
          nsatdown <- downdata$N2sat_mgL[(1+lag):length(downdata$WaterTemp_C)]
        }
        # For bayesian modeling
        if(modType == "bayes"){
          if(gas == "O2") {
            # perform MCMC
            # see documentation on mcmc
            met.post <- mcmc::metrop(tspost_O2, initial = start, nbatch = nbatch,
                                     scale = scale,tempup = tempup,
                                     tempdown = tempdown, oxyup = oxyup,
                                     osatup = osatup, osatdown = osatdown,
                                     oxydown = oxydown, z = z,
                                     light = light, tt = tt, K600mean = K600mean,
                                     K600sd = K600sd, gas = gas, n = n, debug = TRUE,
                                     nspac = 1, lag = lag)

            # Output plot of random walk
            plot(ts(met.post$batch), main = dateList[i])

            # Calculate overall estimates for each day
            gppr <- quantile(met.post$batch[(2000:nbatch),1], c(0.025, 0.5, 0.975))
            err <- quantile(met.post$batch[(2000:nbatch),2], c(0.025, 0.5, 0.975))
            K600r <- quantile(met.post$batch[(2000:nbatch),3], c(0.025, 0.5, 0.975))
            sr <- quantile(met.post$batch[(2000:nbatch),4], c(0.025, 0.5, 0.975))

            # Add results to vectors
            #date[i] <- ymd(unique(data$date)[i])
            GPP[i] <- gppr[2]
            GPP.lower[i] <- gppr[1]
            GPP.upper[i] <- gppr[3]
            ER[i] <- err[2]
            ER.lower[i] <- err[1]
            ER.upper[i] <- err[3]
            K600[i] <- K600r[2]
            K600.lower[i] <- K600r[1]
            K600.upper[i] <- K600r[3]
            s[i] <- sr[2]
            s.lower[i] <- sr[1]
            s.upper[i] <- sr[3]
            accept[i] <- met.post$accept # log likelihood plus priors, should be about 0.2

            # Create dataframe of predicted metabolism values
            pred.metab <- data.frame(date = dateList, GPP, GPP.lower, GPP.upper, ER, ER.lower,
                                     ER.upper, K600, K600.lower, K600.upper, s, s.lower, s.upper)
          }
          if(gas == "N2"){
            # perform MCMC
            # see documentation on mcmc
            met.post <- mcmc::metrop(tspost_N2, initial = start,
                                     nbatch = nbatch,
                                     scale = scale, tempup = tempup,
                                     tempdown = tempdown, n2up = n2up,
                                     nsatup = nsatup, n2down = n2down,
                                     nsatdown = nsatdown, z = z,
                                     light = light, tt = tt, K600mean = K600mean,
                                     K600sd = K600sd, gas = gas, n = n,
                                     debug = TRUE, nspac = 1,
                                     eqn = eqn, lag = lag)

            if(!grepl(pattern = "blende.+", eqn)){
              # trying to troubleshoot here
              plot(ts(met.post$batch,
                      names = c("NConsume", "DN", "K600", "s")),
                   main = dateList[i])

              # Calculate overall estimates for each day
              nconsumer <- quantile(met.post$batch[(2000:nbatch),1], c(0.025, 0.5, 0.975))
              dnr <- quantile(met.post$batch[(2000:nbatch),2], c(0.025, 0.5, 0.975))
              K600r <- quantile(met.post$batch[(2000:nbatch),3], c(0.025, 0.5, 0.975))
              sr <- quantile(met.post$batch[(2000:nbatch),4], c(0.025, 0.5, 0.975))

              # Add results to vectors
              #date[i] <- ymd(unique(data$date)[i])
              NConsume[i] <- nconsumer[2]
              NConsume.lower[i] <- nconsumer[1]
              NConsume.upper[i] <- nconsumer[3]
              DN[i] <- dnr[2]
              DN.lower[i] <- dnr[1]
              DN.upper[i] <- dnr[3]
              K600[i] <- K600r[2]
              K600.lower[i] <- K600r[1]
              K600.upper[i] <- K600r[3]
              s[i] <- sr[2]
              s.lower[i] <- sr[1]
              s.upper[i] <- sr[3]
              accept[i] <- met.post$accept # log likelihood plus priors, should be about 0.2

              # Create dataframe of predicted metabolism values
              pred.metab <- data.frame(date = dateList, NConsume, NConsume.lower,
                                       NConsume.upper, DN, DN.lower,
                                       DN.upper, K600, K600.lower, K600.upper, s,
                                       s.lower, s.upper)
            }
            if(grepl(pattern = "blended[12]", eqn)){
              plot(ts(met.post$batch,
                      names = c("NOther", "NFix", "DN", #"K600",
                                "s")),
                   main = dateList[i])

              # Calculate overall estimates for each day
              notherr <- quantile(met.post$batch[(2000:nbatch), 1],
                                    c(0.025, 0.5, 0.975))
              nfixr <- quantile(met.post$batch[(2000:nbatch), 2],
                                c(0.025, 0.5, 0.975))
              dnr <- quantile(met.post$batch[(2000:nbatch), 3],
                              c(0.025, 0.5, 0.975))
              #K600r <- quantile(met.post$batch[(2000:nbatch), 4],
               #                 c(0.025, 0.5, 0.975))
              sr <- quantile(met.post$batch[(2000:nbatch), 4],
                             c(0.025, 0.5, 0.975))

              # Add results to vectors
              #date[i] <- ymd(unique(data$date)[i])
              NOther[i] <- notherr[2]
              NOther.lower[i] <- notherr[1]
              NOther.upper[i] <- notherr[3]
              NFix[i] <- nfixr[2]
              NFix.lower[i] <- nfixr[1]
              NFix.upper[i] <- nfixr[3]
              DN[i] <- dnr[2]
              DN.lower[i] <- dnr[1]
              DN.upper[i] <- dnr[3]
              #K600[i] <- K600r[2]
              #K600.lower[i] <- K600r[1]
              #K600.upper[i] <- K600r[3]
              s[i] <- sr[2]
              s.lower[i] <- sr[1]
              s.upper[i] <- sr[3]
              accept[i] <- met.post$accept # log likelihood plus priors, should be about 0.2

              # Create dataframe of predicted metabolism values
              pred.metab <- data.frame(date = dateList, NOther, NOther.lower,
                                       NOther.upper, NFix, NFix.lower, NFix.upper,
                                       DN, DN.lower,
                                       DN.upper, #K600, K600.lower, K600.upper,
                                       s,
                                       s.lower, s.upper)
            } # Close blended
            if(eqn == "blended3"){
              plot(ts(met.post$batch,
                      names = c("NOther", "NFix", "DN", "K600",
                                "s")),
                   main = dateList[i])

              # Calculate overall estimates for each day
              notherr <- quantile(met.post$batch[(2000:nbatch), 1],
                                  c(0.025, 0.5, 0.975))
              nfixr <- quantile(met.post$batch[(2000:nbatch), 2],
                                c(0.025, 0.5, 0.975))
              dnr <- quantile(met.post$batch[(2000:nbatch), 3],
                              c(0.025, 0.5, 0.975))
              K600r <- quantile(met.post$batch[(2000:nbatch), 4],
                               c(0.025, 0.5, 0.975))
              sr <- quantile(met.post$batch[(2000:nbatch), 4],
                             c(0.025, 0.5, 0.975))

              # Add results to vectors
              #date[i] <- ymd(unique(data$date)[i])
              NOther[i] <- notherr[2]
              NOther.lower[i] <- notherr[1]
              NOther.upper[i] <- notherr[3]
              NFix[i] <- nfixr[2]
              NFix.lower[i] <- nfixr[1]
              NFix.upper[i] <- nfixr[3]
              DN[i] <- dnr[2]
              DN.lower[i] <- dnr[1]
              DN.upper[i] <- dnr[3]
              K600[i] <- K600r[2]
              K600.lower[i] <- K600r[1]
              K600.upper[i] <- K600r[3]
              s[i] <- sr[2]
              s.lower[i] <- sr[1]
              s.upper[i] <- sr[3]
              accept[i] <- met.post$accept # log likelihood plus priors, should be about 0.2

              # Create dataframe of predicted metabolism values
              pred.metab <- data.frame(date = dateList, NOther, NOther.lower,
                                       NOther.upper, NFix, NFix.lower, NFix.upper,
                                       DN, DN.lower,
                                       DN.upper, K600, K600.lower, K600.upper,
                                       s,
                                       s.lower, s.upper)
            } # Close blended3
          } # Close N2
          } # close if modType == bayes
        } # close else length data > lag
  } # Close looping through each day of datelist NOTE: this may be an unnecessary loop, as only 1 day of data is presented to twostationpostsum at a time

  # Call O2TimeSeries function to return modeled O2 values based on median GPP and ER
  # modeling results
  if(gas == "O2"){
    modeledGas <- O2TimeSeries(GPP = pred.metab$GPP, ER = pred.metab$ER,
                               timeup = timeup, timedown = timedown,
                               oxyup = oxyup, z = z,
                               light = light, tt = tt,
                               tempup = tempup, K600mean = K600mean,
                               gas = gas, n = n, osatup = osatup,
                               osatdown = osatdown, lag = lag,
                               oxydown = oxydown)
  }
  if(gas == "N2"){
    if(!grepl(pattern = "blende.+", eqn)){
      modeledGas <- N2TimeSeries(NConsume = pred.metab$NConsume, DN = pred.metab$DN,
                                 timeup = timeup, timedown = timedown,
                                 n2up = n2up, z = z, light = light, tt = tt,
                                 tempup = tempup, K600mean = K600mean, gas = gas,
                                 n = n, nsatup = nsatup, nsatdown = nsatdown,
                                 lag = lag, n2down = n2down, eqn = eqn)
    }
    if(grepl(pattern = "blende.+", eqn)){
      modeledGas <- N2TimeSeries(NOther = pred.metab$NOther,
                                 NFix = pred.metab$NFix,
                                 DN = pred.metab$DN,
                                 timeup = timeup, timedown = timedown,
                                 n2up = n2up, z = z, light = light, tt = tt,
                                 tempup = tempup, K600mean = K600mean, gas = gas,
                                 n = n, nsatup = nsatup, nsatdown = nsatdown,
                                 lag = lag, n2down = n2down, eqn = eqn)
    }

  }

  # Create output list to return to user
  output <- list(pred.metab = pred.metab, accept = accept, modeledGas = modeledGas)

  return(output)
}
