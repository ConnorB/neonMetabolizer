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
  rawData <- rawData %>%
    dplyr::filter(modelingStrategy == "twoStation")
  # Create 15 min sequence from start to end time of series
  sequence <- data.frame(DateTime_UTC = seq(min(rawData$DateTime_UTC),
                                            max(rawData$DateTime_UTC),
                                            by = "15 min"))
  rawData <- dplyr::full_join(rawData, sequence, by = "DateTime_UTC")
  # Arrange by date
  rawData <- dplyr::arrange(rawData, DateTime_UTC)
  # Split time into chunks of all non-NA DO data or all NA DO data
  modelBlocks <- split(rawData, cumsum(c(TRUE, diff(is.na(rawData$DO_mgL)) != 0)))

  # if nested list has no entries, remove from list
  #if modeling block contains NAs in DO column (therefore, is a section of NA data)
  for (i in seq_along(modelBlocks)) {
    if (any(is.na(modelBlocks[[i]]$DO_mgL))) {
      modelBlocks[[i]] <- NULL
    }
  } #this works properly but throws a subscript out of bounds error, ive played with using appply functions insteadt but i can't get them to work

  #### Pass each modeling block to modeling functions #########################
  # Initiate empty lists to store all the modeling results in
  i = 1
  results <- list()
  results_accept <- list()
  results_modeledO2 <- list()

  for (i in seq_along(modelBlocks)){
    # Initiate empty lists to store each day's worth of modeling results
    model_accept <- list()
    model_metab <- list()
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
                                     Ksd = K600_sd+15,
                                     upName = "S1", downName = "S2",
                                     nbatch = 1e4, scale = 0.3)

      #### Add day of modeled data to list ####################################
      model_metab[[j]] <- metab_out$pred.metab
      model_accept[[j]] <- metab_out$accept
      model_O2[[j]] <- metab_out$oxymodel

      message(modelDate," modeled.")
    } # Exit loop for modeling each day within the modelBlock

    #### Add
    # Combine results for each modelBlock
    results[[i]] <- model_metab %>% dplyr::bind_rows()
    results_accept[[i]] <- unlist(model_accept)
    results_modeledO2[[i]] <- model_O2 %>% dplyr::bind_rows()

  } # Exit loop for modeling the modelBlock

  # Combine metabolism results dataframe
  results <- results %>% dplyr::bind_rows()
  results$accept <- unlist(results_accept)
  results_modeledO2 <- results_modeledO2 %>% dplyr::bind_rows()

  # Return to user
  return(list(results = results,
              modeledO2 = results_modeledO2,
              rawData = rawData))
}

