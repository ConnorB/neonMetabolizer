#' Clean NEON metabolism parameter data and prep for two-station stream metabolism modeling
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
clean_NEON <-function(data){
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
    dataframe$horizontalPosition <- unique(na.omit(dataframe$horizontalPosition))
    #dataframe$referenceLongitude <- unique(na.omit(dataframe$referenceLongitude))
    return(dataframe)
  }
  # Run function
  rawData_S1 <- equalbreaks(rawData_S1)
  rawData_S2 <- equalbreaks(rawData_S2)
  
  #### Fill gaps in raw data ##################################################
  # List of numeric columns where missing data may be present
  numCols <- c("specificConductance_uScm", "DO_mgL", "seaLevelDissolvedOxygenSat",
               "localDissolvedOxygenSat", 
               "pH",
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
  
  #### Return to user #########################################################
  return(rawData)
}

