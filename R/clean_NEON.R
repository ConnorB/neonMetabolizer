#' Clean NEON metabolism parameter data and prep for two-station stream metabolism modeling
#'
#' @importFrom magrittr %>%
#'
#' @param data Dataframe of metabolism parameter time series, `data` in list output of the @request_NEON function
#' @param k600_clean Dataframe of summarized K600 estimates and parameters used for K600 calculations, `k600_clean` in list output of the @request_NEON function
#' @param k600_fit Linear model output for the relationship between mean Q and K600 at the site, `k600_fit` in list output of the @request_NEON function
#'
#' @return
#'
#' @seealso
#'
#' @examples
#'
#' @export
clean_NEON <-function(data, k600_clean, k600_fit){
  # Arrange by date
  data <- dplyr::arrange(data, DateTime_UTC)

  #### Convert from UTC time to solar time #####################################
  # Convert from UTC to solar time
  data$solarTime <-
    streamMetabolizer::convert_UTC_to_solartime(data$DateTime_UTC,
                                                longitude = data$referenceLongitude)

  #### Convert solar time to chron object ######################################
  # Convert solarTime column to chron object
  date <- lubridate::date(data$solarTime)
  time <- paste(lubridate::hour(data$solarTime),
                lubridate::minute(data$solarTime),
                lubridate::second(data$solarTime), sep = ":")
  message("Solar time calculated from local time.")

  # Convert to chron object
  data$dtime <- chron::chron(dates = as.character(date),
                             times = as.character(time),
                             format = c(dates = "y-m-d", times = "h:m:s"))

  #### Search for obviously erronious sensor data ##############################
  # List of sensor data columns
  sensCols <- c("specificConductance_uScm", "DO_mgL", "pH","chlorophyll_ugL",
                "turbidity_NTU", "Nitrate_uMolL", "WaterTemp_C", "AirPres_kPa",
                "WaterPres_kPa", "Depth_m", "Discharge_m3s","Light_PAR")

  # PAR: all nighttime PAR that's below 0 should be set to == 0.
  # Sometimes, esp. if large gaps in data/start of dataset, ARIMA will predict
  # PAR data as dipping below 0, but this is not possible, so hard set to
  # 0 if negative
  # If numeric data (which should all be > 0, is < 0, replace with NA) &
  # if there are any infinite values (Inf), replace these with NA
  data <-
    data %>%
    mutate_at(sensCols, ~ifelse(. < 0, NA, .)) %>%
    mutate_at(sensCols, ~ifelse(is.infinite(.), NA, .))
  message("Obviously erroneous sensor data (readings < 0) eliminated.")

  #### Create equal time breaks from start to end of data series ##############
  # Split data into two dataframes broken up by station
  rawData_S1 <- data %>% dplyr::filter(horizontalPosition == "S1")
  rawData_S2 <- data %>% dplyr::filter(horizontalPosition == "S2")

  # Create vector containing all 15 min chunks from start to end time of series
  sequence <- data.frame(DateTime_UTC = seq(min(data$DateTime_UTC, na.rm = T),
                                            max(data$DateTime_UTC, na.rm = T),
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
    return(dataframe)
  }
  # Run function
  rawData_S1 <- equalbreaks(rawData_S1)
  rawData_S2 <- equalbreaks(rawData_S2)

  #### Fill NA values ##########################################################
  # Fill NAs in numeric columns in each dataframe
  imputeNA_NEON <- function(dataframe){
    # Impute missing values in numeric columns
    for (i in seq_along(sensCols)){
      # Impute gaps using ARIMA model
      imp <- imputeTS::na_kalman(dataframe %>% dplyr::select(sensCols[i]),
                                 model = "auto.arima")
      # Set values in dataframe to imputed values
      dataframe[sensCols[i]] <- imp
    }
    return(dataframe)
  }

  # Run the downstream and upstream sites through NA imputation
  # ARIMA is kinda slow, this might take a couple minutes
  message("Beginning imputation of missing sensor data at site S1. \n   This may take a few minutes.")
  rawData_S1 <- imputeNA_NEON(rawData_S1)
  message("Missing sensor data from site S1 imputed via ARIMA model.")
  message("Beginning imputation of missing sensor data at site S2. \n   This may take a few minutes.")
  rawData_S2 <- imputeNA_NEON(rawData_S2)
  message("Missing sensor data from site S2 imputed via ARIMA model.")
  # Unite data back into a single dataframe
  data <- dplyr::full_join(rawData_S1, rawData_S2)

  #### Add K based on lm relationship #########################################
  predVar <- data.frame(meanQ = data$Discharge_m3s)
  predK600 <- predict.lm(k600_fit, newdata = predVar, interval = "prediction")
  data$K600 <- predK600[,"fit"]
  # Convert from 95% confidence interval to SD
  data$K600_sd <- sqrt(length(k600_clean$k600)) *
    (predK600[,"upr"] - predK600[,"lwr"]) / 3.92
  message("K600 calculated for each sensor timestep based on k600_fit relationship.")

  #### Use Garcia&Gordon and Benson&Krause to calculate DO saturation % ########
  # Calculating DO % saturation under local conditions using equations from
  # the following sources:
  # Garcia, H., and L. Gordon (1992), \emph{Oxygen solubility in seawater:
  # Better fitting equations}, Limnology and Oceanography, 37(6).
  #
  # Benson, B. B. & Krause, D. (1984). \emph{The concentration and isotopic
  # fractionation of oxygen dissolved in freshwater and seawater in equilibrium
  # with the atmosphere.} Limnology and Oceanography, 29(3), 620-632.
  # doi:10.4319/lo.1984.29.3.0620
  #
  # (FYI: This is the same approach as used within the mimsy package)
  #
  # Vapor pressure correction: use the Antoine equation to calculate vapor
  # pressure of water [bar] See NIST Chemistry WebBook for general tables,
  # these parameters valid for temperatures between -18 to 100C (Stull 1947)
  vapor.press <- exp(4.6543 - (1435.264/((data$WaterTemp_C + 273.15) + -64.848)))
  vapor.press <- vapor.press * 0.98692  # conversion from [bar] to [atm]
  # Convert pressure (kPa) to pressure in atm
  AirPres_atm <- data$AirPres_kPa / 101
  # pressure correction [atm] = (current pressure - vapor pressure) /
  # (standard pressure [atm] - vapor pressure)
  press.corr <- (AirPres_atm - vapor.press)/(1 - vapor.press)
  # O2 saturation calculation Combined fit coefficients [umol/kg]
  # (Garcia and Gordon 1992, Table 1)
  A0 <- 5.80818
  A1 <- 3.20684
  A2 <- 4.1189
  A3 <- 4.93845
  A4 <- 1.01567
  A5 <- 1.41575
  B0 <- -7.01211 * 10^-3
  B1 <- -7.25958 * 10^-3
  B2 <- -7.93334 * 10^-3
  B3 <- -5.54491 * 10^-3
  C0 <- -1.32412 * 10^-7
  # Scaled temperature (Garcia and Gordon 1992, eqn. 8)
  TS <- log((298.15 - data$WaterTemp_C)/(273.15 + data$WaterTemp_C))  # log() == natural log (ln)
  # Salinity [per mille], all NEON aquatic sites are freshwater, so set = 0
  S <- 0
  # Calculate O2 saturation concentration at temperature and salinity
  # (Garcia and Gordon 1992, eqn. 8)
  lnO2.sat <- A0 + A1 * TS + A2 * TS^2 + A3 * TS^2 + A3 * TS^3 + A4 *
    TS^4 + A5 * TS^5 + S * (B0 + B1 * TS + B2 * TS^2 + B3 * TS^3) + C0 * S^2
  O2.sat <- exp(lnO2.sat)
  # Correct O2 saturation with pressure correction, solubility.conc units
  # [umol/kg]
  DOsat_uMol.kg <- O2.sat * press.corr
  # Unit conversion: Convert from uMol kg-1 to mg L-1
  # 1 kg H2O = 1 L H2O
  data$DOsat_mgL <- DOsat_uMol.kg * 10^(-6) * 32 * 10^3
  # Calculate percent saturation
  data$DOsat_perc <- (data$DO_mgL / data$DOsat_mgL) * 100
  # Remove localDissolvedOxygenSat column
  data <- select(data, -c(seaLevelDissolvedOxygenSat, localDissolvedOxygenSat))
  message("O2 saturation (mg/L and %) calculated from DO, air pressure, and temperature \n   using Garcia & Gordon (1992) and Benson & Krause (1984) solubility equations.")

  #### Keep only 15-minute dates ###############################################
  # Some sensor measurements are taken on a ~30 sec delay (ex. 12:15:30 at S2
  # and 12:15:00 at S1)
  # Keep only the data where the seconds are 0
  data <- data[second(data$DateTime_UTC) == 0,]

  #### Organize dataframe ######################################################
  data <- data %>%
    relocate(c(solarTime, dtime), .after = c(DateTime_UTC)) %>%
    relocate(c(DOsat_mgL, DOsat_perc), .after = DO_mgL) %>%
    relocate(referenceLongitude, .after = horizontalPosition)

  #### Get travel time between stations #######################################
  #traveltime_fit <- lm(travelTime ~ meanQ, data = k600_clean)

  #data$travelTime_s <- predict.lm(traveltime_fit, newdata = predVar) #CHECK travel time reported in seconds

  # If k600 is outside the bounds of tracer measurements made at the site (flow rate & k600 is higher than any conditions when measurement was taken)
  # travel time might be reported as negative, as the best fit line dips below y = 0 axis. Therefore,
  # If travel time is a negative number, change to 10 seconds - a fast but biologically possible travel time
  #data <- data %>% dplyr::mutate(travelTime_s = ifelse(travelTime_s <= 0, 10, travelTime_s))

  #### Return to user #########################################################
  return(data)
}

