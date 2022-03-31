#' tspost FUNCTION: Calculate posterior probablility of model
#' This function calculates the posterior probability of the model given
#' parameters. Internal function, called within twostationpostsum
#'
#' ARGUMENTS:
#' MET          Dataframe name of cleaned raw two station data (ex. "TS_S1S2")
#' tempup        Temperature data from upstream station, deg C
#' tempdown      Temperature data from downstream station, deg C
#' oxyup         Oxygen data from upstream station, mg-O2/L
#' oxydown       Oxygen data from downstream station, mg-O2/L
#' Light         any light unit
#' tt            travel time, days
#' Kmean
#' Ksd
tspost <- function(MET, tempup, tempdown, oxyup, oxydown, light, tt, z, osat,
                   Kmean, Ksd, ...){
  # Assign the paramters we solve for to easy to understand values
  GPP <- MET[1]
  ER <- MET[2]
  K <- MET[3]
  # Always model the log of variance so that one does not get a
  # negative standard deviation
  sigma <- exp(MET[4])

  lag <- as.numeric(round(tt/0.0104166667))

  metab <- vector(mode = "numeric", length = length(oxyup)) #create empty vector

  # Below is equation 4 in the paper, solving for downstream O2 at each
  # interval. It references other functions:  Kcor converts K600 to KO2
  # for a given temperature.
  for (i in 1:length(oxyup)){
    metab[i] <- (oxyup[i] + ((GPP/z)*(sum(light[i:(i+lag)])/sum(light))) +
                   ER*tt/z +
                   (Kcor(tempup[i],Kmean))*tt*(osat[i] -
                                                 oxyup[i] +
                                                 osat[i])/2) /
      (1+ Kcor(tempup[i],Kmean)*tt/2)
  }

  # likelhood is below.  dnorm caculates the probablity density of a normal
  # distribution, note log.
  loglik <- sum(dnorm(oxydown, metab, sigma, log=TRUE))

  # Priors, note wide distributions for GPP and ER
  prior <- (dnorm(GPP, mean=3.1, sd=6.0, log=TRUE)) +
    (dnorm(ER, mean=-7.1, sd=7.1, log=TRUE)) +
    #(dnorm(K, mean=Kmean, sd=Ksd, log=TRUE)) +
    (dnorm(K, mean=Kmean, sd=Ksd, log=TRUE))

  return(loglik + prior)
}
