#' Model two-station stream metabolism
#'
#' @importFrom magrittr %>%
#'
#' @param data The output from @request_NEON function, a list of 4 objects: `data`: dataframe containing raw site data, `k600_clean`: Dataframe of summarized K600 estimates and parameters used for calculations,`k600_fit`: Linear model output for the relationship between mean Q and K600 at the site, `k600_expanded`: dataframe of all data used in K600 calculations
#' @param nbatch Numeric, number of batches for Bayesian chains. Default is 100,000.
#'
#' @return
#'
#' @seealso
#'
#' @examples
#'
#' @export
fit_twostation <-function(data, nbatch = 1e5){
  # Parse out imput data
  rawData <- data$data
  k600_data <- data$k600_clean
  k600_fit <- data$k600_fit

  #### Add K based on lm relationship #########################################
  predVar <- data.frame(meanQ = rawData$Discharge_m3s)
  rawData$K600 <- predict.lm(k600_fit, newdata = predVar)

  #### Get travel time between stations #######################################
  traveltime_fit <- lm(travelTime ~ meanQ, data = k600_data)
  rawData$travelTime_s <- predict.lm(traveltime_fit, newdata = predVar) #CHECK travel time reported in seconds

  #### Subset data where two station modeling is possible #####################
  # IN FUTURE: experiment with filling short gaps in data
  rawData <- rawData %>%
    dplyr::filter(modelingStrategy == "twoStation")
  # Create vector containing all 15 min chunks from start to end time of series
  sequence <- data.frame(DateTime_UTC = seq(min(rawData$DateTime_UTC),
                                            max(rawData$DateTime_UTC),
                                            by = "15 min"))
  rawData <- dplyr::full_join(rawData, sequence, by = "DateTime_UTC")
  # Arrange by date
  rawData <- dplyr::arrange(rawData, DateTime_UTC)
  # Split time into chunks of all non-NA DO data or all NA DO data
  modelBlocks <- split(rawData, cumsum(c(TRUE, diff(is.na(rawData$DO_mgL)) != 0)))

  # If model block has no entries, remove from list
  for (i in seq_along(modelBlocks)) {
    if (any(is.na(modelBlocks[[i]]$DO_mgL))) {
      modelBlocks[[i]] <- NULL
    }
  } #this works correctly but throws a subscript out of bounds error, ive played with using appply functions insteadt but i can't get them to work

  #### Pass each modeling block to modeling functions #########################
  # Initiate empty lists to store all the modeling results in
  i = 1
  results_accept <- list(rep(NA, length(modelBlocks)))
  results_metab_date <- list(rep(NA, length(modelBlocks)))
  results_metab_GPP <- list(rep(NA, length(modelBlocks)))
  results_metab_GPP.lower <- list(rep(NA, length(modelBlocks)))
  results_metab_GPP.upper <- list(rep(NA, length(modelBlocks)))
  results_metab_ER <- list(rep(NA, length(modelBlocks)))
  results_metab_ER.lower <- list(rep(NA, length(modelBlocks)))
  results_metab_ER.upper <- list(rep(NA, length(modelBlocks)))
  results_metab_K <- list(rep(NA, length(modelBlocks)))
  results_metab_K.lower <- list(rep(NA, length(modelBlocks)))
  results_metab_K.upper <- list(rep(NA, length(modelBlocks)))
  results_metab_s <- list(rep(NA, length(modelBlocks)))
  results_metab_s.lower <- list(rep(NA, length(modelBlocks)))
  results_metab_s.upper <- list(rep(NA, length(modelBlocks)))
  results_modeledO2 <- list()

  for (i in seq_along(modelBlocks)){
    # Initiate empty lists to store each day's worth of modeling results
    # Pre-allocating size of the lists will save computation time, and
    # allocation time for data.frames in R is slow so you should avoid these
    # in loops
    model_accept <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_date <- list(rep(NA, length(unique(lubridate::date(modelBlocks[[i]]$solarTime)))))
    model_metab_GPP <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_GPP.lower <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_GPP.upper <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_ER <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_ER.lower <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_ER.upper <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_K <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_K.lower <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_K.upper <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_s <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_s.lower <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    model_metab_s.upper <- vector(mode = "numeric", length = length(unique(lubridate::date(modelBlocks[[i]]$solarTime))))
    ## Note..... need to figure out how to pre-allocate length of model_02
    model_O2 <- list()

    for (j in seq_along(unique(lubridate::date(modelBlocks[[i]]$solarTime)))){ # For each day within the modelBlock
      # Grab date
      modelDate <- unique(lubridate::date(modelBlocks[[i]]$solarTime))[j]

      # Subset data for selected date
      data_subset <- modelBlocks[[i]][lubridate::date(modelBlocks[[i]]$solarTime) == modelDate,]

      # If not a full day of data, skip
      # Readings every 15 min = 96 readings / day * 2 stations = 192 rows
      if(sum(lubridate::date(data_subset$solarTime) == modelDate) < 192) {
        message("NOTE: Date ", modelDate," contains less than a full day of obervations and was skipped")
        #### Add day of modeled data to list ####################################
        model_accept[[j]] <- NA
        model_metab_date[[j]] <- modelDate
        model_metab_GPP[[j]] <- NA
        model_metab_GPP.lower[[j]] <- NA
        model_metab_GPP.upper[[j]] <- NA
        model_metab_ER[[j]] <- NA
        model_metab_ER.lower[[j]] <- NA
        model_metab_ER.upper[[j]] <- NA
        model_metab_K[[j]] <- NA
        model_metab_K.lower[[j]] <- NA
        model_metab_K.upper[[j]] <- NA
        model_metab_s[[j]] <- NA
        model_metab_s.lower[[j]] <- NA
        model_metab_s.upper[[j]] <- NA
        next
      }
      #### Get mean travel time on selected day, convert from seconds to days #
      travelTime_days <- mean(data_subset$travelTime_s) / 86400

      #### Get mean depth on selected day #####################################
      depth_m <- mean(data_subset$Depth_m)

      #### Get mean and standard deviation of K600 on selected day ############
      K600_mean <- mean(data_subset$K600)
      K600_sd <- sd(data_subset$K600)

      #### Run 2-station metabolism modeling function for selected day ########
      # "Start" denotes the initial state of the markov chain for GPP, ER,
      # K, and s, 'start' is not the same an informative prior
      metab_out <- twostationpostsum(start = c(10, -10, 7, -2.2),
                                     O2data = data_subset,
                                     z = depth_m,
                                     tt = travelTime_days,
                                     Kmean = K600_mean,
                                     Ksd = K600_sd,
                                     upName = "S1", downName = "S2",
                                     nbatch = 1e4, scale = 0.3)

      #### Add day of modeled data to list ####################################
      model_accept[[j]] <- metab_out$accept
      model_metab_date[[j]] <- metab_out$pred.metab$date
      model_metab_GPP[[j]] <- metab_out$pred.metab$GPP
      model_metab_GPP.lower[[j]] <- metab_out$pred.metab$GPP.lower
      model_metab_GPP.upper[[j]] <- metab_out$pred.metab$GPP.upper
      model_metab_ER[[j]] <-  metab_out$pred.metab$ER
      model_metab_ER.lower[[j]] <- metab_out$pred.metab$ER.lower
      model_metab_ER.upper[[j]] <- metab_out$pred.metab$ER.upper
      model_metab_K[[j]] <- metab_out$pred.metab$K
      model_metab_K.lower[[j]] <- metab_out$pred.metab$K.lower
      model_metab_K.upper[[j]] <- metab_out$pred.metab$K.upper
      model_metab_s[[j]] <- metab_out$pred.metab$s
      model_metab_s.lower[[j]] <- metab_out$pred.metab$s.lower
      model_metab_s.upper[[j]] <- metab_out$pred.metab$s.upper

      model_O2[[j]] <- metab_out$oxymodel

      message(modelDate," modeled.")
    } # Exit loop for modeling each day within the modelBlock

    #### Add
    # Combine results for each modelBlock
    results_accept[[i]] <- model_accept
    results_metab_date[[i]] <- model_metab_date
    results_metab_GPP[[i]] <- model_metab_GPP
    results_metab_GPP.lower[[i]] <- model_metab_GPP.lower
    results_metab_GPP.upper[[i]] <- model_metab_GPP.upper
    results_metab_ER[[i]] <- model_metab_ER
    results_metab_ER.lower[[i]] <- model_metab_ER.lower
    results_metab_ER.upper[[i]] <- model_metab_ER.upper
    results_metab_K[[i]] <- model_metab_K
    results_metab_K.lower[[i]] <- model_metab_K.lower
    results_metab_K.upper[[i]] <- model_metab_K.upper
    results_metab_s[[i]] <- model_metab_s
    results_metab_s.lower[[i]] <- model_metab_s.lower
    results_metab_s.upper[[i]] <- model_metab_s.upper
    results_modeledO2[[i]] <- model_O2
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
                        K = unlist(results_metab_K),
                        K.lower = unlist(results_metab_K.lower),
                        K.upper = unlist(results_metab_K.upper),
                        s = unlist(results_metab_s),
                        s.lower = unlist(results_metab_s.lower),
                        s.upper = unlist(results_metab_s.upper))
  results$accept <- unlist(results_accept)
  results_modeledO2 <- dplyr::bind_rows(results_modeledO2)

  # Return to user
  return(list(results = results,
              modeledO2 = results_modeledO2,
              rawData = rawData))
}

