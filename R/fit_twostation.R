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

  #### Create equal time breaks from start to end of data series ##############
  # Split data into two dataframes broken up by station
  rawData_S1 <- rawData %>% dplyr::filter(horizontalPosition == "S1")
  rawData_S2 <- rawData %>% dplyr::filter(horizontalPosition == "S2")

  # Create vector containing all 15 min chunks from start to end time of series
  sequence <- data.frame(DateTime_UTC = seq(min(rawData$DateTime_UTC, na.rm = T),
                                            max(rawData$DateTime_UTC, na.rm = T),
                                            by = "15 min"))
  # Function to merge equal breaks with dataframes
  equalbreaks <- function(dataframe){
    # Join sequence of datetimes with raw data at each station
    dataframe <- dplyr::full_join(dataframe, sequence, by = "DateTime_UTC")
    # Arrange by date
    dataframe <- dplyr::arrange(dataframe, DateTime_UTC)
    # Fill in missing site names and locations from adding missing datetimes
    dataframe$siteID <- unique(na.omit(dataframe$siteID))
    dataframe$horizontalPosition <-unique(na.omit(dataframe$horizontalPosition))
    dataframe$referenceLongitude <-unique(na.omit(dataframe$referenceLongitude))
    return(dataframe)
  }
  # Run function
  rawData_S1 <- equalbreaks(rawData_S1)
  rawData_S2 <- equalbreaks(rawData_S2)

  #### Fill gaps in raw data ##################################################
  # List of numeric columns where missing data may be present
  numCols <- c("specificConductance_uScm", "DO_mgL", "DOsat_pct", "pH",
               "chlorophyll_ugL", "turbidity_NTU", "Nitrate_uMolL", "WaterTemp_C",
               "AirPres_kPa", "WaterPres_kPa", "Depth_m", "Discharge_m3s",
               "Light_PAR")

  # Internal function to fill NAs in numeric columns in each dataframe
  imputeNA_NEON <- function(dataframe){
    # If there is any missing value within row, add a warning that this value
    # has been imputed
    dataframe$warning <-
      ifelse(rowSums(is.na(dataframe %>% dplyr::select(all_of(numCols)))) > 0,
             "This date had one or more missing data points. These were interpolated using an ARIMA model.",
             NA)
    # Impute missing values in numeric columns
    for (i in seq_along(numCols)){
      # Impute gaps using ARIMA model
      imp <- imputeTS::na_kalman(dataframe %>% dplyr::select(numCols[i]),
                                 model = "auto.arima")
      # Set values in dataframe to imputed values
      dataframe[numCols[i]] <- imp
    }
    return(dataframe)
  }

  # Run the downstream and upstream sites through NA imputation
  # ARIMA is kinda slow, this might take a couple minutes
  rawData_S1 <- imputeNA_NEON(rawData_S1)
  message("Missing sensor data from site S1 imputed via ARIMA model.")
  rawData_S2 <- imputeNA_NEON(rawData_S2)
  message("Missing sensor data from site S2 imputed via ARIMA model.")

  # Unite rawData back into a single dataframe
  rawData <- dplyr::full_join(rawData_S1, rawData_S2)
  # Arrange by date
  rawData <- dplyr::arrange(rawData, DateTime_UTC)

  #### Add solar time and convert to chron object #############################
  # Convert from UTC to solar time
  rawData$solarTime <-
    streamMetabolizer::convert_UTC_to_solartime(rawData$DateTime_UTC,
                                                longitude = rawData$referenceLongitude)
  # Convert solarTime column to chron object
  rawData <- rawData %>% dplyr::mutate(date = lubridate::date(solarTime),
                                       time = paste(lubridate::hour(solarTime),
                                              lubridate::minute(solarTime),
                                              lubridate::second(solarTime), sep = ":"))
  # Convert to chron object
  rawData$dtime <- chron::chron(dates = as.character(rawData$date),
                                times = as.character(rawData$time),
                                format = c(dates = "y-m-d", times = "h:m:s"))

  #### Add K based on lm relationship #########################################
  predVar <- data.frame(meanQ = rawData$Discharge_m3s)
  predK600 <- predict.lm(k600_fit, newdata = predVar, interval = "prediction")
  rawData$K600 <- predK600[,"fit"]
  # Convert from 95% confidence interval to SD
  rawData$K600_sd <- sqrt(length(k600_data$k600)) *
    (predK600[,"upr"] - predK600[,"lwr"]) / 3.92

  #### Get travel time between stations #######################################
  traveltime_fit <- lm(travelTime ~ meanQ, data = k600_data)

  rawData$travelTime_s <- predict.lm(traveltime_fit, newdata = predVar) #CHECK travel time reported in seconds

  # If k600 is outside the bounds of tracer measurements made at the site (flow rate & k600 is higher than any conditions when measurement was taken)
  # travel time might be reported as negative, as the best fit line dips below y = 0 axis. Therefore,
  # If travel time is a negative number, change to 10 seconds - a fast but biologically possible travel time
  rawData <- rawData %>% dplyr::mutate(travelTime_s = ifelse(travelTime_s <= 0, 10, travelTime_s))

  #### Subset data where two station modeling is possible #####################
  # IN FUTURE: experiment with filling short gaps in data

  # Split time into chunks of all non-NA DO data or all NA DO data
  #modelBlocks <- split(rawData, cumsum(c(TRUE, diff(is.na(rawData$DO_mgL)) != 0)))

  # If model block has no entries, remove from list
 # for (i in seq_along(modelBlocks)) {
  #  if (any(is.na(modelBlocks[[i]]$DO_mgL))) {
   #   modelBlocks[[i]] <- NULL
    #}
  #} #this works correctly but throws a subscript out of bounds error, ive played with using appply functions insteadt but i can't get them to work

  #### Pass each modeling block to modeling functions #########################
  # Initiate empty lists to store all the modeling results in
  i = 1
  dateList <- unique(rawData$date)
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
    data_subset <- dplyr::filter(rawData, date == modelDate)

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
    K600_mean <- mean(data_subset$K600)
    K600_sd <- mean(data_subset$K600_sd)

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
                        K = unlist(results_metab_K),
                        K.lower = unlist(results_metab_K.lower),
                        K.upper = unlist(results_metab_K.upper),
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
              rawData = rawData))
}

