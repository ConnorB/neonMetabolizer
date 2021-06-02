#' Request data products from NEON API for two-station metabolism modeling
#'
#' \code{request_NEON} communicates with NEON API to request DO,
#'   water temperature, PAR, discharge, depth, and water quality data for the
#'   time period and monitoring station(s) of interest (using NEONScience's
#'   \code{neonUtilities::loadByProduct}). DO
#'   percent saturation is corrected for local elevation
#'   using NEONScience's \code{localPressureDO::calcBK_eq},
#'   and local time is converted to solar time using
#'   \code{streamMetabolizer::convert_UTC_to_solartime}).
#'
#'   \code{request_NEON} additionally requests salt and conductivity slug data
#'   necessary for reaeration rate (K) calculations, and displays an
#'   interactive GUI to the user to identify peaks in slug transport using
#'   NEONScience's \code{reaRate::def.calc.reaeration}. Afterwards,
#'   K600 values are calculated based on the user's input.
#'
#'   Returned is a list containing:
#'   (1) \code{data}: a formatted dataframe with all raw data necessary for
#'   the user to model two-station stream metabolism,
#'   (2) \code{k600_clean}:describe here,
#'   (3) \code{k600_fit}: describe here,
#'   (4) \code{k600_expanded}: describe here.
#'
#' @importFrom magrittr %>%
#'
#' @param NEONsites character string specifying 4-letter NEON site code to
#'   request data from (ex. `"HOPB"`). Can be more than one site
#'   (ex.`c("HOPB", "BLDE")` but be warned data pull will take longer)
#' @param startdate YYYY-MM character string defining start year and monthfor
#'    data request
#' @param enddate YYYY-MM character string defining end year and month for
#'    data request
#'
#' @return List of four dataframes, `data`, `k600_clean`, `k600_fit`, and
#'    `k600_expanded`. NOTE: Explain in detail what these dataframes contain.
#'
#' @seealso \url{https://github.com/NEONScience/NEON-utilities/neonUtilities/}
#'    for details on \code{neonUtilities} package,
#'    \url{https://github.com/NEONScience/NEON-water-quality/localPressureDO}
#'    for details on \code{localPressureDO} package,
#'    \url{https://github.com/NEONScience/NEON-reaeration/reaRate} for details
#'    on \code{reaRate} package.
#'
#' @examples
#' \dontrun{
#' data <- request_NEON(NEONsites = "HOPB", startdate = "2018-01", enddate = "2018-12")
#' }
#'
#' @export
request_NEON <- function(NEONsites, startdate, enddate){
  #### Input parameters ######################################################
  # Define parameters of interest necessary for metabolism modeling
  params <- c("DP1.20288.001", "DP1.20053.001", "DP1.00024.001",
              "DP1.20033.001", "DP4.00130.001")
  names(params) <- c("WaterQual", "Temp", "PAR",
                     "NO3", "Discharge")

  #### Pull NEON data from api ##################################################
  for (i in seq_along(params)){
    # Set dpID to i-th parameter
    dpID <- params[i]

    # Pull NEON data for parameter, saving the data to the variable name from
    # names(params)
    # Pull basic only to save time
    assign(names(dpID),
           value = neonUtilities::loadByProduct(dpID = dpID, site = NEONsites,
                                                startdate = startdate, enddate = enddate,
                                                package = "basic", check.size = F))
  }

  # Pull reaeration data from full period of record
  Reaeration <- neonUtilities::loadByProduct(dpID = "DP1.20190.001",
                                             site = NEONsites,
                                             package = "expanded",
                                             check.size = F)
  FieldDischarge <- neonUtilities::loadByProduct(dpID = "DP1.20048.001",
                                                 site = NEONsites,
                                                 package = "expanded",
                                                 check.size = F)

  #### Correct dissolved oxygen percent saturation ##############################
  # Pull DO-related products for sites of interest
  DOdata <- localPressureDO::getAndFormatData(siteName = NEONsites,
                                              startDate = startdate,
                                              endDate = enddate)
  # Calculate DO percent saturation at local using the Benson-Krause equation,
  # Which is the same method as used by the USGS since 2011
  DOcalcd <- localPressureDO::calcBK_eq(DOdata)
  # Trim DOcalcd dataset
  DOcalcd <-
    DOcalcd %>%
    dplyr::select(horizontalPosition, startDateTime, siteID,
                  dissolvedOxygenSatCorrected)
  # Remove DO percent saturation from WaterQual
  # Merge corrected DO percent saturation with WaterQual dataset
  WaterQual$waq_instantaneous <-
    dplyr::right_join(WaterQual$waq_instantaneous, DOcalcd,
                      by = c("siteID", "startDateTime", "horizontalPosition"))

  #### Format and merge dataframes ##############################################
  # S1 = upstream, S2 = downstream
  ### Nitrate ###
  NO3_data <-
    NO3$NSW_15_minute %>%
    dplyr::select(siteID, startDateTime, surfWaterNitrateMean) %>%
    dplyr::mutate(startDateTime = lubridate::with_tz(startDateTime,
                                                     tz = "UTC")) %>%
    dplyr::rename(DateTime_UTC = startDateTime,
                  Nitrate_uMolL = surfWaterNitrateMean)

  ### Water quality ###
  WQ_data <-
    WaterQual$waq_instantaneous %>%
    dplyr::select(siteID, horizontalPosition, startDateTime,
                  specificConductance, dissolvedOxygen,
                  dissolvedOxygenSatCorrected, pH, chlorophyll, turbidity) %>%
    dplyr::mutate(horizontalPosition = dplyr::recode(horizontalPosition,
                                                     "101" = "S1",
                                                     "102" = "S2"),
                  startDateTime = lubridate::with_tz(startDateTime,
                                                     tz = "UTC")) %>%
    dplyr::rename(DateTime_UTC = startDateTime, DO_mgL = dissolvedOxygen,
                  DOsat_pct = dissolvedOxygenSatCorrected,
                  chlorophyll_ugL = chlorophyll, turbidity_NTU = turbidity,
                  specificConductance_uScm = specificConductance) %>%
    dplyr::filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))

  ### Temp ###
  Temp_data <-
    Temp$TSW_5min %>%
    dplyr::select(siteID, horizontalPosition, startDateTime,
                  surfWaterTempMean) %>%
    dplyr::mutate(horizontalPosition = dplyr::recode(horizontalPosition,
                                                     "101" = "S1",
                                                     "102" = "S2"),
                  startDateTime = lubridate::with_tz(startDateTime,
                                                     tz = "UTC")) %>%
    dplyr::rename(DateTime_UTC = startDateTime,
                  WaterTemp_C = surfWaterTempMean) %>%
    dplyr::filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))

  ### Air pressure ###
  AirPres_data <-
    BP_1min %>%
    dplyr::select(siteID, startDateTime, staPresMean) %>%
    dplyr::mutate(startDateTime = lubridate::with_tz(startDateTime,
                                                     tz = "UTC")) %>%
    dplyr::rename(DateTime_UTC = startDateTime, AirPres_kPa = staPresMean) %>%
    dplyr::filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))

  ### Discharge & Depth ###
  Discharge_data <-
    Discharge$csd_continuousDischarge %>%
    dplyr::select(siteID, endDate, calibratedPressure, equivalentStage,
                  maxpostDischarge) %>%
    dplyr::mutate(endDate = lubridate::with_tz(endDate, tz = "UTC"),
                  maxpostDischarge = maxpostDischarge / 1000) %>% # convert discharge, measured as L/s, to m3/s
    dplyr::rename(DateTime_UTC = endDate, WaterPres_kPa = calibratedPressure,
                  Depth_m = equivalentStage, Discharge_m3s = maxpostDischarge) %>%
    dplyr::filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))

  ### Light ###
  PAR_data <-
    PAR$PARPAR_1min %>%
    dplyr::select(siteID, startDateTime, PARMean) %>%
    dplyr::mutate(startDateTime = lubridate::with_tz(startDateTime,
                                                     tz = "UTC")) %>%
    dplyr::rename(DateTime_UTC = startDateTime, Light_PAR = PARMean) %>%
    dplyr::filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))

  ### Merge all into data dataframe ###
  data <- dplyr::full_join(WQ_data, NO3_data, by = c("siteID","DateTime_UTC"))
  data <- dplyr::left_join(data, Temp_data,
                           by = c("siteID","horizontalPosition","DateTime_UTC"))
  data <- dplyr::left_join(data, AirPres_data, by = c("siteID","DateTime_UTC"))
  data <- dplyr::left_join(data, Discharge_data, by = c("siteID","DateTime_UTC"))
  data <- dplyr::left_join(data, PAR_data, by = c("siteID","DateTime_UTC"))

  ### Convert datetime to chron solartime ###
  # Grab longitude of currently in use sensor stations
  sensorPos <-
    WaterQual$sensor_positions_20288 %>%
    dplyr::filter(end == "") %>%
    dplyr::select(siteID, HOR.VER, referenceLongitude) %>%
    dplyr::mutate(HOR.VER = regmatches(HOR.VER, regexpr(pattern = "^\\d{3}",
                                                        text = HOR.VER)),
                  HOR.VER = dplyr::recode(HOR.VER, "101" = "S1", "102" = "S2")) %>%
    dplyr::rename(horizontalPosition = HOR.VER)

  # Add longitude for use in solartime conversion
  data <- dplyr::right_join(data, sensorPos, by = c("siteID","horizontalPosition"))

  # Convert from UTC to solar time
  data$solarTime <- streamMetabolizer::convert_UTC_to_solartime(data$DateTime_UTC,
                                                                longitude = data$referenceLongitude)

  #################### Reaeration Rate (K) Calculations #########################
  # Format reaeration data product
  Reaeration_data <-
    reaRate::def.format.reaeration(rea_backgroundFieldCondData =
                                     Reaeration$rea_backgroundFieldCondData,
                                   rea_backgroundFieldSaltData =
                                     Reaeration$rea_backgroundFieldSaltData,
                                   rea_fieldData = Reaeration$rea_fieldData,
                                   rea_plateauMeasurementFieldData =
                                     Reaeration$rea_plateauMeasurementFieldData,
                                   rea_externalLabDataSalt =
                                     Reaeration$rea_externalLabDataSalt,
                                   rea_externalLabDataGas =
                                     Reaeration$rea_externalLabDataGas,
                                   rea_widthFieldData =
                                     Reaeration$rea_widthFieldData,
                                   dsc_fieldData = FieldDischarge$dsc_fieldData,
                                   dsc_individualFieldData =
                                     FieldDischarge$dsc_individualFieldData)

  k600_expanded <-
    reaRate::def.calc.reaeration(inputFile = Reaeration_data,
                               loggerData = Reaeration$rea_conductivityFieldData,
                               namedLocation = "namedLocation",
                               injectionTypeName = "injectionType",
                               eventID = "eventID",
                               stationToInjectionDistance = "stationToInjectionDistance",
                               plateauGasConc = "plateauGasConc",
                               corrPlatSaltConc = "corrPlatSaltConc",
                               hoboSampleID = "hoboSampleID",
                               discharge = "fieldDischarge",
                               waterTemp = "waterTemp",
                               wettedWidth = "wettedWidth",
                               plot = TRUE,
                               savePlotPath = NULL,
                               processingInfo = NULL)

  k600 <- k600_expanded$outputDF
  k600_clean <- k600[k600$k600 > 0 & k600$travelTime > 0,]

  # Linear model of Q vs K600
  lmk600 <- lm(k600 ~ meanQ, data = k600_clean)

  ################### Output data to user #######################################
  output <- list(data = data,
                 k600_clean = k600_clean,
                 k600_fit = lmk600,
                 k600_expanded = k600_expanded)

  return(output)
}
