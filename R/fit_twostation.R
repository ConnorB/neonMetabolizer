#' Model two-station stream metabolism
#'
#' @importFrom magrittr %>%
#'
#' @param data The output from @request_NEON function, a list of 4 objects: `data`: dataframe containing raw site data, `k600_clean`: Dataframe of summarized K600 estimates and parameters used for calculations,`k600_fit`: Linear model output for the relationship between mean Q and K600 at the site, `k600_expanded`: dataframe of all data used in K600 calculations
#' @param modType Model type as a character string. Either "mle" for maximum likelyhood estimation or "bayes" for Bayesian estimation.
#' @param nbatch Numeric, number of batches for Bayesian chains. Default is 100,000.
#'
#' @return
#'
#' @seealso
#'
#' @examples
#'
#' @export
fit_twostation <-function(data, modType, nbatch = 1e5, gas, n = 0.5,
                          K600_mean = NULL, K600_sd = NULL,
                          upstreamName = "S1", downstreamName = "S2"){
  # Add date column to dataframe based on solar time - this gets passed into
  # twostationpostsum
  data$date <- lubridate::date(data$solarTime)

  #### Subset data where two station modeling is possible #####################
  # IN FUTURE: experiment with filling short gaps in data

  # Split time into chunks of all non-NA DO data or all NA DO data
  #modelBlocks <- split(data, cumsum(c(TRUE, diff(is.na(data$DO_mgL)) != 0)))

  # If model block has no entries, remove from list
 # for (i in seq_along(modelBlocks)) {
  #  if (any(is.na(modelBlocks[[i]]$DO_mgL))) {
   #   modelBlocks[[i]] <- NULL
    #}
  #} #this works correctly but throws a subscript out of bounds error, ive played with using appply functions insteadt but i can't get them to work

  #### Pass each modeling block to modeling functions #########################
  # Initiate empty lists to store all the modeling results in
  i = 1
  dateList <- unique(data$date)
  results_accept <- list(rep(NA, length(dateList)))
  results_metab_date <- list(rep(NA, length(dateList)))
  results_metab_K600 <- list(rep(NA, length(dateList)))
  results_metab_K600.lower <- list(rep(NA, length(dateList)))
  results_metab_K600.upper <- list(rep(NA, length(dateList)))
  results_metab_s <- list(rep(NA, length(dateList)))
  results_metab_s.lower <- list(rep(NA, length(dateList)))
  results_metab_s.upper <- list(rep(NA, length(dateList)))
  results_warnings <- list(rep(NA, length(dateList)))
  results_modeledGas <- list()

  if(gas == "O2"){
    results_metab_GPP <- list(rep(NA, length(dateList)))
    results_metab_GPP.lower <- list(rep(NA, length(dateList)))
    results_metab_GPP.upper <- list(rep(NA, length(dateList)))
    results_metab_ER <- list(rep(NA, length(dateList)))
    results_metab_ER.lower <- list(rep(NA, length(dateList)))
    results_metab_ER.upper <- list(rep(NA, length(dateList)))
  }
  if(gas == "N2"){
    results_metab_NConsume <- list(rep(NA, length(dateList)))
    results_metab_NConsume.lower <- list(rep(NA, length(dateList)))
    results_metab_NConsume.upper <- list(rep(NA, length(dateList)))
    results_metab_DN <- list(rep(NA, length(dateList)))
    results_metab_DN.lower <- list(rep(NA, length(dateList)))
    results_metab_DN.upper <- list(rep(NA, length(dateList)))
  }
  i <- 1
  for (i in seq_along(dateList)){
    # Grab date
    modelDate <- dateList[i]
    # Subset data for selected date
    data_subset <- dplyr::filter(data, date == modelDate)

    # If there is a discharge 0 reading during day, skip day of data, tell user
    # stream is dry or intermittently dry on this day
    if (any(data_subset$Discharge_m3s == 0)) {
      message("NOTE: Stream is dry or intermittently dry (discharge = 0 m3s) on ", modelDate, ". This day was skipped.")
      # Add day of modeled data to list with NAs
      results_accept[i] <- NA
      results_metab_date[i] <- modelDate
      results_metab_GPP[i] <- NA
      results_metab_GPP.lower[i] <- NA
      results_metab_GPP.upper[i] <- NA
      results_metab_ER[i] <- NA
      results_metab_ER.lower[i] <- NA
      results_metab_ER.upper[i] <- NA
      results_metab_K600[i] <- NA
      results_metab_K600.lower[i] <- NA
      results_metab_K600.upper[i] <- NA
      results_metab_s[i] <- NA
      results_metab_s.lower[i] <- NA
      results_metab_s.upper[i] <- NA
      results_warnings[i] <- "Stream dry or intermittently dry on date"
      next
    }

    #### Get mean travel time on selected day, convert from seconds to days #
    travelTime_days <- mean(data_subset$travelTime_s) / 86400

    #### Get mean depth on selected day #####################################
    depth_m <- mean(data_subset$Depth_m)

    #### Reaeration ###########################################################
    # Get mean and standard deviation of K600 on selected day
    if(is.null(K600_mean) & is.null(K600_sd)){
      K600_mean <- mean(data_subset$K600)
      K600_sd <- mean(data_subset$K600_sd)
    }

    #### Run 2-station metabolism modeling function for selected day ########
    # "Start" denotes the initial state of the markov chain for GPP, ER,
    # K, and s, 'start' is not the same an informative prior
    metab_out <- twostationpostsum(start = c(3.1, -7.1, 7, -2.2),
                                   modType = modType,
                                   data = data_subset,
                                   z = depth_m,
                                   tt = travelTime_days,
                                   K600mean = K600_mean,
                                   K600sd = K600_sd,
                                   upName = upstreamName,
                                   downName = downstreamName,
                                   gas = gas, n = n,
                                   nbatch = nbatch, scale = 0.3)

    #### Add day of modeled data to list ####################################
    results_accept[i] <- metab_out$accept
    results_metab_date[i] <- metab_out$pred.metab$date
    results_metab_K600[i] <- metab_out$pred.metab$K600
    results_metab_K600.lower[i] <- metab_out$pred.metab$K600.lower
    results_metab_K600.upper[i] <- metab_out$pred.metab$K600.upper
    results_metab_s[i] <- metab_out$pred.metab$s
    results_metab_s.lower[i] <- metab_out$pred.metab$s.lower
    results_metab_s.upper[i] <- metab_out$pred.metab$s.upper
    results_warnings[i] <- NA

    if(gas == "N2"){
      results_metab_NConsume[i] <- metab_out$pred.metab$NConsume
      results_metab_NConsume.lower[i] <- metab_out$pred.metab$NConsume.lower
      results_metab_NConsume.upper[i] <- metab_out$pred.metab$NConsume.upper
      results_metab_DN[i] <-  metab_out$pred.metab$DN
      results_metab_DN.lower[i] <- metab_out$pred.metab$DN.lower
      results_metab_DN.upper[i] <- metab_out$pred.metab$DN.upper
    }
    if(gas == "O2"){
      results_metab_GPP[i] <- metab_out$pred.metab$GPP
      results_metab_GPP.lower[i] <- metab_out$pred.metab$GPP.lower
      results_metab_GPP.upper[i] <- metab_out$pred.metab$GPP.upper
      results_metab_ER[i] <-  metab_out$pred.metab$ER
      results_metab_ER.lower[i] <- metab_out$pred.metab$ER.lower
      results_metab_ER.upper[i] <- metab_out$pred.metab$ER.upper
    }

    results_modeledGas[[i]] <- metab_out$modeledGas

    message(modelDate," modeled.")
  } # Exit loop for modeling the modelBlock

  if(gas == "N2"){
    results <- data.frame(date = as.Date(unlist(results_metab_date),
                                         origin = "1970-01-01"),
                          NConsume = unlist(results_metab_NConsume),
                          NConsume.lower = unlist(results_metab_NConsume.lower),
                          NConsume.upper = unlist(results_metab_NConsume.upper),
                          DN = unlist(results_metab_DN),
                          DN.lower = unlist(results_metab_DN.lower),
                          DN.upper = unlist(results_metab_DN.upper),
                          K600 = unlist(results_metab_K600),
                          K600.lower = unlist(results_metab_K600.lower),
                          K600.upper = unlist(results_metab_K600.upper),
                          s = unlist(results_metab_s),
                          s.lower = unlist(results_metab_s.lower),
                          s.upper = unlist(results_metab_s.upper),
                          warnings = unlist(results_warnings)
    )
  }
  if(gas == "O2"){
    # Combine metabolism results dataframe
    results <- data.frame(date = as.Date(unlist(results_metab_date),
                                         origin = "1970-01-01"),
                          GPP = unlist(results_metab_GPP),
                          GPP.lower = unlist(results_metab_GPP.lower),
                          GPP.upper = unlist(results_metab_GPP.upper),
                          ER = unlist(results_metab_ER),
                          ER.lower = unlist(results_metab_ER.lower),
                          ER.upper = unlist(results_metab_ER.upper),
                          K600 = unlist(results_metab_K600),
                          K600.lower = unlist(results_metab_K600.lower),
                          K600.upper = unlist(results_metab_K600.upper),
                          s = unlist(results_metab_s),
                          s.lower = unlist(results_metab_s.lower),
                          s.upper = unlist(results_metab_s.upper),
                          warnings = unlist(results_warnings)
    )
  }

  results$accept <- unlist(results_accept)
  results_modeledGas <- dplyr::bind_rows(results_modeledGas)

  # Return to user
  return(list(results = results,
              modeledGas = results_modeledGas,
              data = data))
}

