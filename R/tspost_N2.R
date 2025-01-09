#' \code{tspost} Calculate posterior probablility of two-station metabolism model.
#'
#' This function calculates the posterior probability of the the two station metabolism model given
#' parameters. Internal function, called within twostationpostsum
#'
#' @param MET Dataframe name of cleaned raw two station data (ex. "TS_S1S2")
#' @param tempup Temperature data from upstream station, deg C
#' @param tempdown Temperature data from downstream station, deg C
#' @param n2up Oxygen data from upstream station, mg-O2/L
#' @param n2down Oxygen data from downstream station, mg-O2/L
#' @param light any light unit
#' @param tt travel time, days
#' @param z Description
#' @param nsatup Description
#' @param nsatdown Description
#' @param Kmean Description
#' @param Ksd Description
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
tspost_N2 <- function(MET, tempup, tempdown, n2up, n2down, light, tt, z, nsatup,
                      nsatdown, K600mean, K600sd, gas, n){
  # Assign the parameters we solve for to easy to understand values
  NConsume <- MET[1]
  DN <- MET[2]
  K600 <- MET[3]
  # Always model the log of variance so that one does not get a
  # negative standard deviation
  sigma <- exp(MET[4])

  lag <- as.numeric(round(tt/0.0104166667))
  # If lag is really small (less than the time between sample collections), set lag = 1
  if(lag < 1){
    lag <- 1
  }

  metab <- vector(mode = "numeric", length = length(n2up)) #create empty vector

  # Solving for downstream N2 at each interval. It references other functions:
  # Kcor converts K600 to KO2 for a given temperature & gas type.
  for (i in 1:length(n2up)){
    # (this function from nifong et al )
    metab[i] <-
      (n2up[i] + (
        (NConsume/z) *
          ( sum(light[i:(i+lag)]) / sum(light) )
        ) + DN * tt/z +
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

  # likelhood is below.  dnorm calculates the probablity density of a normal
  # distribution, note log.
  loglik <- sum(dnorm(n2down, metab, sigma, log=TRUE))

  prior <-
    # Priors for NConsume and DN based on Nifong et al 2020
    (dnorm(NConsume, mean = -0.1, sd = 5, log=TRUE)) +
    (dnorm(DN, mean = 0.1, sd = 5, log=TRUE)) +
    # Priors for NConsume and DN based on median N fix in Marcarelli et al 2024
    # and median Denitrification rates in Marcarelli et al 2024 & LINXII
   # (dnorm(NConsume, mean = -4.104199e-05, sd = 0.002207237, log=TRUE)) +
    #(dnorm(DN, mean = 0.01640903, sd = 0.0866865, log=TRUE)) +
    #(dnorm(K, mean=Kmean, sd=Ksd, log=TRUE)) +
    (dnorm(K600, mean=K600mean, sd=K600sd, log=TRUE))

  return(loglik + prior)
}
