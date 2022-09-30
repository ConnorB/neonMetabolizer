#' \code{twostationpostsum} Calculate two station MCMC or Bayes during one day.
#'
#' Internal function. Function runs the MCMC and returns posterior distributions
#'
#' @param O2data Dataframe containing formatted raw two station data, as returned by request_NEON in $data
#' @param upName Character string denoting name of upstream station (ex. "S1")
#' @param downName Character string denoting name of downstream station (ex. "S2")
#' @param start Description
#' @param z Numeric, average depth on day of modeling
#' @param tt Numeric, travel time between stations on day of modeling
#' @param Kmean Mean K for informative prior
#' @param Ksd Standard deviation of K for informative prior
#' @param nbatch Number of MCMC trials
#' @param scale Description
#' @param modType Description
#'
#' @returns Populate here
#'
#' @references Populate here
#'
#' @example Populate here
twostationpostsum <- function(O2data, upName, downName, start, z, tt, Kmean, Ksd,
                              nbatch, scale, modType) {
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
        osat <- updata$DOsat_mgL[1:as.numeric(length(updata$WaterTemp_C)-lag)]
        tempdown <- downdata$WaterTemp_C[(1+lag):length(downdata$WaterTemp_C)]
        oxydown <- downdata$DO_mgL[(1+lag):length(downdata$WaterTemp_C)]
        light <- downdata$Light_PAR

        # For MLE modeling
        if(modType == "mle"){
          # Wrap within trycatch, which will print a descriptive error if nlm fails
          tryCatch(
            {
              # Try to evaluate nlm with data
              met.post <- nlm(tspost,
                              hessian = TRUE, control=list(trace=TRUE, maxit=2000),
                              p = start, tempup = tempup, tempdown = tempdown,
                              oxyup = oxyup, osat = osat, oxydown = oxydown, z = z,
                              light = light, tt = tt, Kmean = Kmean, Ksd = Ksd)

              # Compute standard errors
              nlmErr <- sqrt(diag(solve(-met.post[["hessian"]])))
                # If this can't compute - it will return NaN

              # 95% confidence intervals will be parameter +- 1.96*std error
              GPP[i] <- met.post$estimate[1]
              GPP.lower[i] <- met.post$estimate[1] - 1.96*nlmErr[1]
              GPP.upper[i] <- met.post$estimate[1] + 1.96*nlmErr[1]
              ER[i] <- met.post$estimate[2]
              ER.lower[i] <- met.post$estimate[2] - 1.96*nlmErr[2]
              ER.upper[i] <- met.post$estimate[2] + 1.96*nlmErr[2]
              K[i] <- met.post$estimate[3]
              K.lower[i] <- met.post$estimate[3] - 1.96*nlmErr[3]
              K.upper[i] <- met.post$estimate[3] + 1.96*nlmErr[3]
              s[i] <- met.post$estimate[4]
              s.lower[i] <- met.post$estimate[4] - 1.96*nlmErr[4]
              s.upper[i] <- met.post$estimate[4] - 1.96*nlmErr[4]
              accept[i] <- met.post$code # See nlm documentation for code information
            },
            error = function(err){
              # If nlm returns an error, print this error message to user
              message(paste("ERROR: Non-finite estimate supplied by `nlm`. Modeling on date",
                      dateList[i], "failed."))

              # Add NA values to dataframe
              GPP[i] <- NA
              GPP.lower[i] <- NA
              GPP.upper[i] <- NA
              ER[i] <- NA
              ER.lower[i] <- NA
              ER.upper[i] <- NA
              K[i] <- NA
              K.lower[i] <- NA
              K.upper[i] <- NA
              s[i] <- NA
              s.lower[i] <- NA
              s.upper[i] <- NA
              accept[i] <- "Model failure, non-finite estimate supplied by `nlm`."
            }, # Close error argument
            warning = function(warn){
              # If nlm returns a warning, print this error message to user
              message(paste0("WARNING: NA/Inf replaced by maximum positive value \n\tduring one or more timesteps on ",
                             dateList[i], "."))
            }, #Close warning argument
            finally = {
              # Add data to dataframe of predicted metabolism values
              pred.metab <- data.frame(date = dateList, GPP, GPP.lower,
                                       GPP.upper, ER, ER.lower, ER.upper, K,
                                       K.lower, K.upper, s, s.lower, s.upper,
                                       accept)
            } # Close finally argument
          ) # Close tryCatch
        }
        # For bayesian modeling
        if(modType == "bayes"){
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

          # Create dataframe of predicted metabolism values
          pred.metab <- data.frame(date = dateList, GPP, GPP.lower, GPP.upper, ER, ER.lower,
                                   ER.upper, K, K.lower, K.upper, s, s.lower, s.upper)
          } # close if modType == bayes
        } # close else length data > lag
  } # Close looping through each day of datelist NOTE: this may be an unnecessary loop, as only 1 day of data is presented to twostationpostsum at a time

  # Call O2TimeSeries function to return modeled O2 values based on median GPP and ER
  # modeling results
  oxymodel <- O2TimeSeries(GPP = pred.metab$GPP, ER = pred.metab$ER,
                           O2data = O2data, Kmean = Kmean, z = z, tt = tt,
                           upName = upName, downName = downName)

  # Create output list to return to user
  output <- list(pred.metab = pred.metab, accept = accept, oxymodel = oxymodel)

  return(output)
}
