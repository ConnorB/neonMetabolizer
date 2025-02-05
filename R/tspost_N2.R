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
                      nsatdown, K600mean, K600sd, gas, n, eqn, lag){
  # Assign the parameters we solve for to easy to understand values
  NConsume <- MET[1]
  DN <- MET[2]
  K600 <- MET[3]
  # Always model the log of variance so that one does not get a
  # negative standard deviation
  sigma <- exp(MET[4])

  metab <- vector(mode = "numeric", length = length(n2up)) #create empty vector

  # Solve for downstream N2 at each interval
  if(eqn == "Nifong_et_al_2020"){
    for (i in 1:length(n2up)){
      # (this function from nifong et al 2020)
      metab[i] <- nifong(i = i, n2up = n2up, NConsume = NConsume, z = z,
                         light = light, DN = DN, tt = tt, tempup = tempup,
                         K600mean = K600mean, gas = gas, n = n, nsatup = nsatup,
                         nsatdown = nsatdown, lag = lag)
    }
  }
  if(eqn == "light_independent"){
    for (i in 1:length(n2up)){
      metab[i] <- lightIndependent(i = i, n2up = n2up, NConsume = NConsume, z = z,
                                   DN = DN, tt = tt, tempup = tempup,
                                   K600mean = K600mean, gas = gas, n = n,
                                   nsatup = nsatup,
                                   nsatdown = nsatdown)
    }
  }
  if(grepl(pattern = "blende.+", eqn)){
    # Assign the parameters we solve for to easy to understand values
    NOther <- MET[1]
    NFix <- MET[2]
    DN <- MET[3]
    #K600 <- MET[4]
    sigma <- exp(MET[4])

    if(eqn == "blended1"){
      for (i in 1:length(n2up)){
        metab[i] <- blended1(i = i, n2up = n2up, NOther = NOther,
                            NFix = NFix, z = z, light = light,
                            DN = DN, tt = tt, tempup = tempup,
                            K600mean = K600mean, gas = gas, n = n,
                            nsatup = nsatup,
                            nsatdown = nsatdown, lag = lag)
      }
    }

    if(eqn == "blended2"){
      for (i in 1:length(n2up)){
        metab[i] <- blended2(i = i, n2up = n2up, NOther = NOther,
                            NFix = NFix, z = z, light = light,
                            DN = DN, tt = tt, tempup = tempup,
                            K600mean = K600mean, gas = gas, n = n,
                            nsatup = nsatup,
                            nsatdown = nsatdown, lag = lag)
      }
    }

      if(eqn == "blended3"){
        K600 <- MET[4]
        sigma <- exp(MET[4])

        for (i in 1:length(n2up)){
          metab[i] <- blended1(i = i, n2up = n2up, NOther = NOther,
                               NFix = NFix, z = z, light = light,
                               DN = DN, tt = tt, tempup = tempup,
                               K600mean = K600mean, gas = gas, n = n,
                               nsatup = nsatup,
                               nsatdown = nsatdown, lag = lag)
        }
      }
  }


  # likelhood is below.  dnorm calculates the probablity density of a normal
  # distribution, note log.
  loglik <- sum(dnorm(n2down, metab, sigma, log=TRUE))

  if(eqn %in% c("Nifong_et_al_2020", "light_independent")){
    prior <-
      # Priors for NConsume and DN based on Kelly lit review
      (dnorm(NConsume, mean = -9.6e-5, sd = 0.30, log=TRUE)) +
      (dnorm(DN, mean = 0.07, sd = 0.25, log=TRUE)) +
      # Priors for NConsume and DN from Nifong et al 2020
      #(dnorm(NConsume, mean = -0.1, sd = 5, log=TRUE)) +
      #(dnorm(DN, mean = 0.1, sd = 5, log=TRUE)) +
      (dlnorm(K600, meanlog=K600mean, sdlog=K600sd, log=TRUE))
  }
  if(grepl(pattern = "blended[12]", eqn)){
    prior <-
      # Priors for NConsume and DN from Nifong et al 2020
      (dnorm(NOther, mean = -0.1, sd = 5, log=TRUE)) +
      (dnorm(NFix, mean = -0.1, sd = 5, log=TRUE)) +
      (dlnorm(DN, meanlog = 0.1, sdlog = 5, log=TRUE)) #+
      #(dlnorm(K600, meanlog=K600mean, sdlog=K600sd, log=TRUE))
  }
  if(eqn == "blended3"){
    prior <-
      # Priors for NConsume and DN from Nifong et al 2020
      (dnorm(NOther, mean = -0.1, sd = 5, log=TRUE)) +
      (dnorm(NFix, mean = -0.1, sd = 5, log=TRUE)) +
      (dlnorm(DN, meanlog = 0.1, sdlog = 5, log=TRUE)) +
      (dlnorm(K600, meanlog=K600mean, sdlog=K600sd, log=TRUE))
  }

  return(loglik + prior)
}
