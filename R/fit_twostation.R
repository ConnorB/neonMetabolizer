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
fit_twostation <-function(data, modType, nbatch = 1e5){
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
  results_metab_GPP <- list(rep(NA, length(dateList)))
  results_metab_GPP.lower <- list(rep(NA, length(dateList)))
  results_metab_GPP.upper <- list(rep(NA, length(dateList)))
  results_metab_ER <- list(rep(NA, length(dateList)))
  results_metab_ER.lower <- list(rep(NA, length(dateList)))
  results_metab_ER.upper <- list(rep(NA, length(dateList)))
  results_metab_K <- list(rep(NA, length(dateList)))
  results_metab_K.lower <- list(rep(NA, length(dateList)))
  results_metab_K.upper <- list(rep(NA, length(dateList)))
  results_metab_s <- list(rep(NA, length(dateList)))
  results_metab_s.lower <- list(rep(NA, length(dateList)))
  results_metab_s.upper <- list(rep(NA, length(dateList)))
  results_warnings <- list(rep(NA, length(dateList)))
  results_modeledO2 <- list()

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
      results_metab_K[i] <- NA
      results_metab_K.lower[i] <- NA
      results_metab_K.upper[i] <- NA
      results_metab_s[i] <- NA
      results_metab_s.lower[i] <- NA
      results_metab_s.upper[i] <- NA
      results_warnings[i] <- "Stream dry or intermittently dry on date"
      next
    }

    # If not a full day of data, skip
    # Readings every 15 min = 96 readings / day * 2 stations = 192 rows
    #if(sum(lubridate::date(data_subset$solarTime) == modelDate) < 192) {
  #    message("NOTE: Date ", modelDate," contains less than a full day of obervations and was skipped")
      #### Add day of modeled data to list ####################################
   #   model_accept[i] <- NA
    #  model_metab_date[i] <- modelDate
     # model_metab_GPP[i] <- NA
      #model_metab_GPP.lower[i] <- NA
      #model_metab_GPP.upper[i] <- NA
      #model_metab_ER[i] <- NA
      #model_metab_ER.lower[i] <- NA
      #model_metab_ER.upper[i] <- NA
      #model_metab_K[i] <- NA
      #model_metab_K.lower[i] <- NA
      #model_metab_K.upper[i] <- NA
      #model_metab_s[i] <- NA
      #model_metab_s.lower[i] <- NA
      #model_metab_s.upper[i] <- NA
      #next
    #}

    #### Get mean travel time on selected day, convert from seconds to days #
    travelTime_days <- mean(data_subset$travelTime_s) / 86400

    #### Get mean depth on selected day #####################################
    depth_m <- mean(data_subset$Depth_m)

    #### Get mean and standard deviation of K600 on selected day ############
    k600_mean <- mean(data_subset$k600)
    k600_sd <- mean(data_subset$k600_sd)
    #k600_sd <- 0.001
    #k600_sd <- mean(data_subset$k600_sd) # will this improve fits? giving 5% wiggle room
    #k600_sd <- 5 # super constrained fit on k

    #### Run 2-station metabolism modeling function for selected day ########
    # "Start" denotes the initial state of the markov chain for GPP, ER,
    # K, and s, 'start' is not the same an informative prior
    metab_out <- twostationpostsum(start = c(3.1, -7.1, 7, -2.2),
                                   modType = modType,
                                   O2data = data_subset,
                                   z = depth_m,
                                   tt = travelTime_days,
                                   Kmean = k600_mean,
                                   Ksd = k600_sd,
                                   upName = "S1", downName = "S2",
                                   nbatch = nbatch, scale = 0.3)

    #### Add day of modeled data to list ####################################
    results_accept[i] <- metab_out$accept
    results_metab_date[i] <- metab_out$pred.metab$date
    results_metab_GPP[i] <- metab_out$pred.metab$GPP
    results_metab_GPP.lower[i] <- metab_out$pred.metab$GPP.lower
    results_metab_GPP.upper[i] <- metab_out$pred.metab$GPP.upper
    results_metab_ER[i] <-  metab_out$pred.metab$ER
    results_metab_ER.lower[i] <- metab_out$pred.metab$ER.lower
    results_metab_ER.upper[i] <- metab_out$pred.metab$ER.upper
    results_metab_K[i] <- metab_out$pred.metab$K
    results_metab_K.lower[i] <- metab_out$pred.metab$K.lower
    results_metab_K.upper[i] <- metab_out$pred.metab$K.upper
    results_metab_s[i] <- metab_out$pred.metab$s
    results_metab_s.lower[i] <- metab_out$pred.metab$s.lower
    results_metab_s.upper[i] <- metab_out$pred.metab$s.upper
    results_warnings[i] <- NA

    results_modeledO2[[i]] <- metab_out$oxymodel

    message(modelDate," modeled.")
  } # Exit loop for modeling the modelBlock

  # Combine metabolism results dataframe
  results <- data.frame(date = as.Date(unlist(results_metab_date),
                                       origin = "1970-01-01"),
                        GPP = unlist(results_metab_GPP),
                        GPP.lower = unlist(results_metab_GPP.lower),
                        GPP.upper = unlist(results_metab_GPP.upper),
                        ER = unlist(results_metab_ER),
                        ER.lower = unlist(results_metab_ER.lower),
                        ER.upper = unlist(results_metab_ER.upper),
                        k = unlist(results_metab_K),
                        k.lower = unlist(results_metab_K.lower),
                        k.upper = unlist(results_metab_K.upper),
                        s = unlist(results_metab_s),
                        s.lower = unlist(results_metab_s.lower),
                        s.upper = unlist(results_metab_s.upper),
                        warnings = unlist(results_warnings)
                        )
  results$accept <- unlist(results_accept)
  results_modeledO2 <- dplyr::bind_rows(results_modeledO2)

  # Return to user
  return(list(results = results,
              modeledO2 = results_modeledO2,
              data = data))
}

