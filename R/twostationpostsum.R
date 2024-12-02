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
                              nbatch, scale, modType, gas, n) {
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
  }

  i <- 1
  # For loop iterating through each day of data
  for (i in 1:length(dateList)){
    # Seperate data into upstream and downstream sections
    updata <- data[data$horizontalPosition == upName,]
    downdata <- data[data$horizontalPosition == downName,]

    # NOTE: In Hall et al, the number below was 0.00347222, which corresponds to:
    # 0.00347222 days = 5 minutes, as their O2 sensors took 5-minute readings.
    # Here, our sensors took 15 minute readings. So:
    # 15 minutes = 0.0104166667 days
    # number of 15 min readings between up and down probe corresponding to travel
    # time tt
    lag <- as.numeric(round(tt/0.0104166667))
    # In the event that travel time is really fast between stations (lag = 0) set lag = 1
    if(lag == 0){
      lag <- 1
    }

    # trim the ends of the oxy and temp data by the lag so that oxydown[1]
    # is the value that is the travel time later than oxy up.  The below calls
    # are designed to work with our data structure.
    if(length(updata$DO_mgL) < lag) {
      # If there is less data during a day than the lag interval, move to next day
      message("ERROR: model not computed for ", dateList[i], " as insufficient observations provided.")
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
      }
      next
      } else{
        # Else continue
        tempup <- updata$WaterTemp_C[1:as.numeric(length(updata$WaterTemp_C)-lag)] # trim the end by the lag
        tempdown <- downdata$WaterTemp_C[(1+lag):length(downdata$WaterTemp_C)]
        lightup <- updata$Light_PAR[1:as.numeric(length(updata$Light_PAR)-lag)]
        lightdown <- downdata$Light_PAR[(1+lag):length(downdata$Light_PAR)]
        lightmean <- (lightup + lightdown)/2
        # Script needs to sum i:i+lag, problem for end of time string.
        finallight <- lightmean[length(lightmean)]
        # Extrapolate final light as last light values
        lightmean <- c(lightmean, rep(finallight, 1+lag))
        light <- lightdown
        if(gas == "O2"){
          oxyup <- updata$DO_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          osatup <- updata$DOsat_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          oxydown <- downdata$DO_mgL[(1+lag):length(downdata$WaterTemp_C)]
          osatdown <- downdata$DOsat_mgL[(1+lag):length(downdata$WaterTemp_C)]
        }
        if(gas == "N2"){
          n2up <- updata$N2_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          nsat <- updata$N2sat_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
          n2down <- downdata$N2_mgL[(1+lag):length(downdata$WaterTemp_C)]
        }
        # For MLE modeling
        if(modType == "mle"){
          # Wrap within trycatch, which will print a descriptive error if nlm fails
          tryCatch(
            {
              # Try to evaluate nlm with data
              if(gas == "O2"){
                met.post <- nlm(tspost,
                                hessian = TRUE, #control=list(trace=TRUE, maxit=2000),
                                gas = "O2", n = 0.5,
                                p = start, tempup = tempup, tempdown = tempdown,
                                oxyup = oxyup, osat = osat, oxydown = oxydown, z = z,
                                light = light, tt = tt, K600mean = K600mean, K600sd = K600sd)
              }
              if(gas == "N2"){
                met.post <- nlm(tspost_N2,
                                hessian = TRUE, control=list(trace=TRUE, maxit=2000),
                                p = start, tempup = tempup, tempdown = tempdown,
                                oxyup = oxyup, osat = osat, oxydown = oxydown, z = z,
                                light = light, tt = tt, K600mean = K600mean, K600sd = K600sd)
              }

              # Compute standard errors
              nlmErr <- sqrt(diag(solve(-met.post[["hessian"]])))
                # If this can't compute - it will return NaN
              # 95% confidence intervals will be parameter +- 1.96*std error
              K600[i] <- met.post$estimate[3]
              K600.lower[i] <- met.post$estimate[3] - 1.96*nlmErr[3]
              K600.upper[i] <- met.post$estimate[3] + 1.96*nlmErr[3]
              s[i] <- met.post$estimate[4]
              s.lower[i] <- met.post$estimate[4] - 1.96*nlmErr[4]
              s.upper[i] <- met.post$estimate[4] - 1.96*nlmErr[4]
              accept[i] <- met.post$code # See nlm documentation for code information

              if(gas == "O2"){
                GPP[i] <- met.post$estimate[1]
                GPP.lower[i] <- met.post$estimate[1] - 1.96*nlmErr[1]
                GPP.upper[i] <- met.post$estimate[1] + 1.96*nlmErr[1]
                ER[i] <- met.post$estimate[2]
                ER.lower[i] <- met.post$estimate[2] - 1.96*nlmErr[2]
                ER.upper[i] <- met.post$estimate[2] + 1.96*nlmErr[2]
              }
              if(gas == "N2"){
                NConsume[i] <- met.post$estimate[1]
                NConsume.lower[i] <- met.post$estimate[1] - 1.96*nlmErr[1]
                NConsume.upper[i] <- met.post$estimate[1] + 1.96*nlmErr[1]
                DN[i] <- met.post$estimate[2]
                DN.lower[i] <- met.post$estimate[2] - 1.96*nlmErr[2]
                DN.upper[i] <- met.post$estimate[2] + 1.96*nlmErr[2]
              }
            },
            error = function(err){
              # If nlm returns an error, print this error message to user
              message(paste("ERROR: Non-finite estimate supplied by `nlm`. Modeling on date",
                      dateList[i], "failed."))
              # Add NA values to dataframe
              K600[i] <- NA
              K600.lower[i] <- NA
              K600.upper[i] <- NA
              s[i] <- NA
              s.lower[i] <- NA
              s.upper[i] <- NA
              accept[i] <- "Model failure, non-finite estimate supplied by `nlm`."
              if(gas == "O2"){
                GPP[i] <- NA
                GPP.lower[i] <- NA
                GPP.upper[i] <- NA
                ER[i] <- NA
                ER.lower[i] <- NA
                ER.upper[i] <- NA
              }
              if(gas == "N2") {
                NConsume[i] <- NA
                NConsume.lower[i] <- NA
                NConsume.upper[i] <- NA
                DN[i] <- NA
                DN.lower[i] <- NA
                DN.upper[i] <- NA
              }
            }, # Close error argument
            warning = function(warn){
              # If nlm returns a warning, print this error message to user
              message(paste0("WARNING: NA/Inf replaced by maximum positive value \n\tduring one or more timesteps on ",
                             dateList[i], "."))
            }, #Close warning argument
            finally = {
              # Add data to dataframe of predicted metabolism values
              if(gas == "O2"){
                pred.metab <- data.frame(date = dateList, GPP, GPP.lower,
                                         GPP.upper, ER, ER.lower, ER.upper, K600,
                                         K600.lower,  K600.upper, s, s.lower, s.upper,
                                         accept)
              } else if(gas == "N2") {
                pred.metab <- data.frame(date = dateList, NConsume, NConsume.lower,
                                         NConsume.upper, DN, DN.lower, DN.upper,K600,
                                         K600.lower,K600.upper, s, s.lower, s.upper,
                                         accept)
              }
            } # Close finally argument
          ) # Close tryCatch
        }
        # For bayesian modeling
        if(modType == "bayes"){
          if(gas == "O2") {
            # perform MCMC
            # see documentation on mcmc
            met.post <- mcmc::metrop(tspost, initial = start, nbatch = nbatch,
                                     scale = scale,tempup = tempup,
                                     tempdown = tempdown, oxyup = oxyup,
                                     osatup = osatup, osatdown = osatdown,
                                     oxydown = oxydown, z = z,
                                     light = light, tt = tt, K600mean = K600mean,
                                     K600sd = K600sd, gas = gas, n = n, debug = TRUE,
                                     nspac = 1)

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
            met.post <- mcmc::metrop(tspost_N2, initial = start, nbatch = nbatch,
                                     scale = scale, tempup = tempup,
                                     tempdown = tempdown, n2up = n2up,
                                     nsat = nsat, n2down = n2down,  z = z,
                                     light = light, tt = tt, K600mean = K600mean,
                                     K600sd = K600sd, gas = gas, n = n, debug = TRUE, nspac = 1)

            # trying to troubleshoot here
            plot(ts(met.post$batch), main = dateList[i])

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
            pred.metab <- data.frame(date = dateList, NConsume, NConsume.lower, NConsume.upper, DN, DN.lower,
                                     DN.upper, K600, K600.lower, K600.upper, s, s.lower, s.upper)
          }

          } # close if modType == bayes
        } # close else length data > lag
  } # Close looping through each day of datelist NOTE: this may be an unnecessary loop, as only 1 day of data is presented to twostationpostsum at a time

  # Call O2TimeSeries function to return modeled O2 values based on median GPP and ER
  # modeling results
  if(gas == "O2"){
    modeledGas <- O2TimeSeries(GPP = pred.metab$GPP, ER = pred.metab$ER,
                               data = data, K600mean = K600mean, z = z, tt = tt,
                               upName = upName, downName = downName, gas = gas,
                               n = n)
  }
  if(gas == "N2"){
    modeledGas <- N2TimeSeries(NConsume = pred.metab$NConsume, DN = pred.metab$DN,
                               data = data, K600mean = K600mean, z = z, tt = tt,
                               upName = upName, downName = downName, gas = gas,
                               n = n)
  }

  # Create output list to return to user
  output <- list(pred.metab = pred.metab, accept = accept, modeledGas = modeledGas)

  return(output)
}
