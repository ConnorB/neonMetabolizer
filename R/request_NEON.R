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
              "DP1.20033.001", "DP4.00130.001", "DP1.20093.001",
              "DP1.20072.001", "DP1.00004.001")
  names(params) <- c("WaterQual", "Temp", "PAR",
                     "NO3", "Discharge", "WaterChem",
                     "AquaticPlants", "Barometer")

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

  # NOTE: as of 18-Jan-2022, Water Quality data product now includes local DO
  # percent saturation. Therefore, correcting dissolved oxygen percent saturation
  # using localPressureDO package is unneeded.

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
                  seaLevelDissolvedOxygenSat,
                  localDissolvedOxygenSat, pH, chlorophyll, turbidity) %>%
    dplyr::mutate(horizontalPosition = dplyr::recode(horizontalPosition,
                                                     "101" = "S1",
                                                     "102" = "S2"),
                  startDateTime = lubridate::with_tz(startDateTime,
                                                     tz = "UTC")) %>%
    dplyr::rename(DateTime_UTC = startDateTime, DO_mgL = dissolvedOxygen,
                  chlorophyll_ugL = chlorophyll, turbidity_NTU = turbidity,
                  specificConductance_uScm = specificConductance) %>%
    dplyr::filter(lubridate::minute(DateTime_UTC) %in% c(0, 15, 30, 45))

  ### Water chem ###
  WC_data <-
    WaterChem$swc_externalLabDataByAnalyte %>%
    dplyr::select(siteID, namedLocation, startDate,
                  analyte, analyteConcentration, analyteUnits) %>%
    dplyr::mutate(namedLocation = stringr::str_extract(string = namedLocation,
                                                       pattern = "\\w\\d$"),
                  analyteUnits = dplyr::recode(analyteUnits,
                                               "milligramsPerLiter" = "_mgL"),
                  analyteUnits = tidyr::replace_na(analyteUnits, ""),
                  analyte = paste0(analyte, analyteUnits)) %>%
    dplyr::select(siteID, namedLocation, startDate, analyte, analyteConcentration) %>%
    tidyr::pivot_wider(names_from = analyte, values_from = analyteConcentration,
                       values_fn = mean) %>%
    dplyr::rename(NO3NO2_mgNL = `NO3+NO2 - N_mgL`, NH4_mgNL = `NH4 - N_mgL`) %>%
    dplyr::mutate(DIN_mgNL = NO3NO2_mgNL + NH4_mgNL,
                  DIN_molNL = DIN_mgNL / 10^6 / 14,
                  TDP_molL = TDP_mgL / 10^6 / 30)

  ### Pull percent cover from Aquatic Plant data ###
  coverData <-
    AquaticPlants$apc_pointTransect %>%
    dplyr::select(siteID, pointNumber, transectDistance, collectDate,
                  habitatType, substrate, remarks) %>%
    dplyr::filter(collectDate >=
                    lubridate::ymd_hms(paste0(startdate, "-01 00:00:00"),
                                       tz = tz(collectDate)) &
                    collectDate <=
                    lubridate::ymd_hms(paste0(enddate, "-01 00:00:00"),
                                       tz = tz(collectDate))) %>%
    dplyr::group_by(siteID, pointNumber) %>%
    dplyr::arrange(transectDistance, .by_group = TRUE)

  percCover_full <-
    coverData %>%
    dplyr::count(substrate) %>%
    dplyr::rename(transect = pointNumber, observationCount = n) %>%
    dplyr::mutate(sumObservations = sum(observationCount),
                  percentCover = observationCount / sumObservations * 100)

  percCover_summary <-
    percCover_full %>%
    dplyr::select(transect, substrate, percentCover) %>%
    tidyr::pivot_wider(names_from = substrate, values_from = percentCover) %>%
    #dplyr::rename(coarseWoodyDebris = CWD, leafLitter = `leaf litter`) %>%
    replace(is.na(.), 0) %>%
    dplyr::ungroup(transect) %>%
    dplyr::select(-transect) %>%
    dplyr::summarise_all(mean)

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
    Barometer$BP_1min %>%
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
  data <- dplyr::full_join(data, Temp_data,
                           by = c("siteID","horizontalPosition","DateTime_UTC"))
  data <- dplyr::full_join(data, AirPres_data, by = c("siteID","DateTime_UTC"))
  data <- dplyr::full_join(data, Discharge_data, by = c("siteID","DateTime_UTC"))
  data <- dplyr::full_join(data, PAR_data, by = c("siteID","DateTime_UTC"))

  ### Convert datetime to chron solartime ###
  # Grab longitude of currently in use sensor stations
  sensorPos <-
    WaterQual$sensor_positions_20288 %>%
    dplyr::filter(referenceEnd == "") %>%
    dplyr::select(siteID, HOR.VER, referenceLongitude) %>%
    dplyr::mutate(HOR.VER = regmatches(HOR.VER, regexpr(pattern = "^\\d{3}",
                                                        text = HOR.VER)),
                  HOR.VER = dplyr::recode(HOR.VER, "101" = "S1", "102" = "S2")) %>%
    dplyr::filter(HOR.VER == "S1" | HOR.VER == "S2") %>%
    dplyr::rename(horizontalPosition = HOR.VER)

  # Add longitude for use in solartime conversion
  data <- dplyr::full_join(data, sensorPos, by = c("siteID","horizontalPosition"))

  #################### Reaeration Rate (K) Calculations #########################
  # Format reaeration data product
  Reaeration_data <-
    reaRate::def.format.reaeration(rea_backgroundFieldCondData =
                                     Reaeration$rea_backgroundFieldCondData,
                                   rea_backgroundFieldSaltData =
                                     Reaeration$rea_backgroundFieldSaltData,
                                   rea_fieldData =
                                     Reaeration$rea_fieldData,
                                   rea_plateauMeasurementFieldData =
                                     Reaeration$rea_plateauMeasurementFieldData,
                                   rea_plateauSampleFieldData =
                                     Reaeration$rea_plateauSampleFieldData,
                                   rea_externalLabDataSalt =
                                     Reaeration$rea_externalLabDataSalt,
                                   rea_externalLabDataGas =
                                     Reaeration$rea_externalLabDataGas,
                                   rea_widthFieldData =
                                     Reaeration$rea_widthFieldData,
                                   dsc_fieldData = FieldDischarge$dsc_fieldData,
                                   dsc_individualFieldData =
                                     FieldDischarge$dsc_individualFieldData,
                                   dsc_fieldDataADCP =
                                     FieldDischarge$dsc_fieldDataADCP)

  # Check if reaeration data is only slug injections
  if(all(Reaeration_data$injectionType %in% c("model","model - slug","model - CRI"))){
    # If every line of the Reaeration_data is a model injection
    message(paste0("All reaeration measurements are model injection, therefore did not calculate gas loss rate."))
    lmk600 <- NA
    k600_clean <- NA
    k600_expanded <- Reaeration_data
  } else {
    # Feed reaeration data to reaeration calculation tool
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
    # Remove K estimates that are less than 0
    k600_clean <- k600[which(k600$k600 > 0),]

    # Linear model of Q vs K600
    lmk600 <- lm(k600 ~ meanQ, data = k600_clean)

    # Plot
    plot(x = k600_clean$meanQ,y = k600_clean$k600,
         main = paste0(unique(k600_clean$siteID),
                       ": Mean Q Against Estimated K600"),
         xlab = "Mean Q", ylab = "k600")
    abline(lmk600, lty = 2, col = "red")
    mtext(paste0("k600 = ", round(summary(lmk600)$coefficients[1], 3), " + ",
                 round(summary(lmk600)$coefficients[2], 3), " * MeanQ, ",
                 "R^2 = ", round(summary(lmk600)$r.squared, 3), side = 3),
          col = "red")
  }

  ################### Output data to user #######################################
  # Remove dataframes from the environment that are generated within the
  # functions
  suppressWarnings(rm(list = c("readme_20288", "sensor_positions_20288",
                               "variables_20288","waq_instantaneous"),
                      envir = .GlobalEnv))

  # Output to user
  output <- list(data = data,
                 waterQual = WC_data,
                 percentCover = percCover_summary,
                 k600_clean = k600_clean,
                 k600_fit = lmk600,
                 k600_expanded = k600_expanded)

  return(output)
}
