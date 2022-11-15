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

  # Save copy of raw data, this will be returned to user in output
  rawData <- data

  #### Convert from UTC time to solar time #####################################
  # Convert from UTC to solar time
  data$solarTime <-
    convert_UTC_to_solartime(data$DateTime_UTC,
                             longitude = data$referenceLongitude)

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
  message("> Obviously erroneous sensor data (readings < 0) eliminated from: \n   ",
          paste(sensCols[1:6], collapse = ", "), ",\n   ",
          paste(sensCols[7:12], collapse = ", "))
  # Check for wonky DO data - a few DO_mgL are in the 600 range, eliminate
  message("> Obviously erroneous DO_mgL (readings > 300) eliminated. This was ",
          sum(na.omit(data$DO_mgL) > 300), " datapoints \n   (or ",
          round(sum(na.omit(data$DO_mgL) > 300)/length(data$DO_mgL)*100, digits = 3),
          "% of datapoints).")
  data <-
    data %>%
    mutate(DO_mgL = case_when(DO_mgL > 300 ~ NA_real_, TRUE ~ DO_mgL))

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
  message("> Beginning imputation of missing sensor data at site S1. \n   This may take a few minutes.")
  rawData_S1 <- imputeNA_NEON(rawData_S1)
  message("> Missing sensor data from site S1 imputed via ARIMA model.")
  message("> Beginning imputation of missing sensor data at site S2. \n   This may take a few minutes.")
  rawData_S2 <- imputeNA_NEON(rawData_S2)
  message("> Missing sensor data from site S2 imputed via ARIMA model.")
  # Unite data back into a single dataframe
  data <- dplyr::full_join(rawData_S1, rawData_S2)

  #### Add K based on lm relationship #########################################
  # Check for discharge == 0
  message("> Discharge (m3 s-1) is equal to 0.00 at ",
          sum(data$Discharge_m3s <= 0), " datapoints \n   (or ",
          round(sum(data$Discharge_m3s <= 0)/length(data$Discharge_m3s)*100, digits = 3),
          "% of datapoints). At these datapoints, \n   discharge has now been set to == 0.0000001 to \n   allow for log transformation.")
  # Change any datapoints where discharge == 0 to discharge == 0.0000001
  # This will allow for happy log transformations at the travel time fit step
  data$Discharge_m3s[data$Discharge_m3s <= 0] <- 0.0000001

  # Check whether k600 data was available from NEON:
  if(is.na(k600_clean)){
    message("> Reaeration measurements (K600) were not available at this site.")
  } else{

    # Get fit statistics for measured K600 vs Q fit
    modSum <- summary(k600_fit)

    # Plot fit relationship for k600
    plot(x = k600_clean$meanQ_cms, y = k600_clean$k600.clean,
         main = "Discharge vs k600", xlab = "Discharge (m3 s-1)",
         ylab = "k600 (d-1)")
    abline(k600_fit, lty = 2, col = "red")
    text(x = max(k600_clean$meanQ_cms)*0.8, y = max(k600_clean$k600.clean)*0.8,
         labels = paste0("R2 = ", format(modSum$adj.r.squared, digits = 3)))
    text(x = max(k600_clean$meanQ_cms)*0.8, y = max(k600_clean$k600.clean)*0.7,
         labels = paste0("p = ", format(modSum$coefficients[2,4], digits = 2)))

    # Check if relationship between K and Q is statistically significant
    # If it is, add K values to dataframe using linear relationship
    # If it's not, use mean of measured K for dataframe K values
    if(modSum$coefficients[2,4] < 0.05){
      message("> Fit of measured k600 vs. Q (p = ", format(modSum$coefficients[2,4], digits = 2),
              ") significant. Therefore, k600 was calculated for each sensor timestep\n   based on k600_fit relationship.")
      # Use measured K600 vs Q fit to predict K600 at each timestep
      predVar <- data.frame(meanQ_cms = data$Discharge_m3s)
      predk600 <- predict.lm(k600_fit, newdata = predVar, interval = "prediction")
      data$k600 <- predk600[,"fit"]
      # Grab 95% confidence interval
      data$k600_lower <- predk600[,"lwr"]
      data$k600_upper <- predk600[,"upr"]
      data$k600_sd <- ((predk600[,"upr"] - predk600[,"lwr"]) / 3.92) *
        sqrt(length(k600_clean$k600.clean)) # Convert from 95% confidence interval to SD
      } else{
        message("> Fit of measured k600 vs. Q (p = ", format(modSum$coefficients[2,4], digits = 2),
                ") insignificant. Therefore, k600 for each sensor timestep\n   set as mean of measured k600 values.")
        # Set k600 to mean of measured k600 values
        data$k600 <- mean(k600_clean$k600.clean, na.rm = TRUE)
        data$k600_sd <- sd(k600_clean$k600.clean, na.rm = TRUE)
      }

    #### Add travel time based on log-lm relationship ############################
    k600_clean$logmeanQ_cms <- log10(k600_clean$meanQ_cms)
    predVar <- data.frame(logmeanQ_cms = log10(data$Discharge_m3s))
    TT_fit <- lm(peakMaxTravelTime ~ logmeanQ_cms, data = k600_clean)
    predTT <- predict.lm(TT_fit, newdata = predVar, interval = "prediction")
    data$travelTime_s <- predTT[,"fit"]
    message("> Travel time between sensor stations calculated for each timestep \n   based on linear relationship between peakMaxTravelTime and Log10(meanQ_cms).")
    # Depending on the fit of the model, there may be some predicted travel times
    # that are less than 0. Obviously this is impossible. To fix this, let's
    # remove less than 0 travel times and replace them with 1 second. This will
    # allow us to see how many travel times have been replaced, also
    message("> ", sum(data$travelTime_s <= 0), " predicted travel times (or ",
            round(sum(data$travelTime_s <= 0)/length(data$travelTime_s)*100,
                  digits = 2),
            "% of datapoints) were predicted as <= 0 \n   according to the log-linear model. As a travel time of <= 0 is \n   impossible, negative travel times have \n   been changed to travel time == 1 second.")
    data$travelTime_s[data$travelTime_s <=0] <- 1 # second
    # Find travel times that are unlikely to be physically possible, such as
    # travel time > 1 day (86400 s). These estimates likely come from periods
    # where the streambed is dry, and flow is super low.
    message("> ", sum(data$travelTime_s >= 86400), " predicted travel times (or ",
            round(sum(data$travelTime_s >= 86400)/length(data$travelTime_s)*100,
                  digits = 2),
            "% of datapoints) were predicted as \n   greater than 1 day. These estimates likely come from periods when the \n   flow is low or stagnant. Travel times exceeding 1 day (86400 sec) \n   have been changed to 0.9 days (77760 sec), so that two-station \n   modeling can be attempted for this day of data.")
    data$travelTime_s[data$travelTime_s >= 86400] <- 77760 # seconds, == 0.9 days
    modSum <- summary(TT_fit)
    # Plot fit relationship for travel time
    plot(x = k600_clean$logmeanQ_cms, y = k600_clean$peakMaxTravelTime,
         main = "Log(Discharge) vs Travel time", xlab = "Log(Discharge (m3 s-1))",
         ylab = "Travel time between stations (s)")
    abline(TT_fit, lty = 2, col = "red")
    text(x = min(k600_clean$logmeanQ_cms)*0.5, y = max(k600_clean$peakMaxTravelTime)*0.8,
         labels = paste0("R2 = ", format(modSum$adj.r.squared, digits = 3)))
    text(x = min(k600_clean$logmeanQ_cms)*0.5, y = max(k600_clean$peakMaxTravelTime)*0.7,
         labels = paste0("p = ", format(modSum$coefficients[2,4], digits = 2)))
  } # Close else k600 available from neon


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
  data$DOsat_pct <- (data$DO_mgL / data$DOsat_mgL) * 100
  # Remove localDissolvedOxygenSat column
  data <- select(data, -c(seaLevelDissolvedOxygenSat, localDissolvedOxygenSat))
  message("> O2 saturation (mg/L and %) calculated from DO, air pressure, and temperature \n   using Garcia & Gordon (1992) and Benson & Krause (1984) solubility equations.")

  #### Keep only 15-minute dates ###############################################
  # Some sensor measurements are taken on a ~30 sec delay (ex. 12:15:30 at S2
  # and 12:15:00 at S1)
  # Keep only the data where the seconds are 0
  data <- data[second(data$DateTime_UTC) == 0,]

  #### Organize cleanData ######################################################
  data <- data %>%
    relocate(c(solarTime), .after = c(DateTime_UTC)) %>%
    relocate(c(DOsat_mgL, DOsat_pct), .after = DO_mgL) %>%
    relocate(referenceLongitude, .after = horizontalPosition)

  #### Format cleanData for metabolism modeling ################################
  # Seperate upstream and downstream data from the dataframe,
  Up <- data %>% filter(horizontalPosition == "S1")
  Down <- data %>% filter(horizontalPosition == "S2")
  # Keep only the data necessary for each sampling site
  # From the upstream dataset, we need Oup and Osatup
  Up <- Up %>%
    distinct(solarTime, .keep_all = TRUE) %>%
    select(DateTime_UTC, solarTime, DO_mgL, DOsat_mgL) %>%
    rename(solar.time = solarTime, DO.obs.up = DO_mgL, DO.sat.up = DOsat_mgL)
  # From the downstream dataset, we need Odown, Osatdown, and all other
  # sensor parameters measured at sensor S2
  Down <- Down %>%
    distinct(solarTime, .keep_all = TRUE) %>%
    select(DateTime_UTC, solarTime, DO_mgL, DOsat_mgL, Depth_m, Discharge_m3s, WaterTemp_C, Light_PAR, k600,
           travelTime_s) %>%
    rename(solar.time = solarTime, DO.obs.down = DO_mgL, DO.sat.down = DOsat_mgL,
           depth = Depth_m, discharge = Discharge_m3s, temp.water = WaterTemp_C,
           light = Light_PAR, tt = travelTime_s)
  # Merge dataframes
  formattedData <- dplyr::full_join(Down, Up, by = "DateTime_UTC") %>%
    rename(solar.time = solar.time.x) %>%
    select(-solar.time.y, -DateTime_UTC) # absolutely no idea why it's not seeing
    # the solar.time columns as equivalent. Super lazy fix for now, come back
    #to this later

  #### Return to user #########################################################
  outList <- list(formattedData = formattedData,
                  cleanData = data,
                  rawData = rawData,
                  k600_clean = k600_clean,
                  k600_fit = k600_fit)

  return(outList)
}
