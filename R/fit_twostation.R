#' Model two-station stream metabolism
#'
#' @importFrom magrittr %>%
#'
#' @param data Dataframe of raw data as returned in the first object from @request_NEON 's output list
#' @param k600_fit Linear model output for the relationship between mean Q and K600 at the site, the third object returned by @request_NEON
#'
#' @return
#'
#' @seealso
#'
#' @examples
#'
#' @export
fit_twostation <-function(data, k600_fit){
  # Parse out imput data
  rawData <- data
  lmK600 <- k600_fit

  #### Add prior for K based on lm results ####################################

  #### Subset data where two station modeling is possible #####################
  rawData <- rawData %>%
    dplyr::filter(modelingStrategy == "twoStation")
  # Create 15 min sequence from start to end time of series
  sequence <- data.frame(DateTime_UTC = seq(min(rawData$DateTime_UTC),
                                            max(rawData$DateTime_UTC),
                                            by = "15 min"))
  rawData <- dplyr::full_join(rawData, sequence, by = "DateTime_UTC")
  # Arrange by date
  rawData <- dplyr::arrange(a, DateTime_UTC)
  # Split time into sections of non NA data
  modelBlocks <- split(a, cumsum(c(TRUE, diff(is.na(a$DO_mgL)) != 0)))
  # if nested list has no entries, remove from list
  for (i in 1:length(modelBlocks)) {
    if (nrow(modelBlocks[[i]]) == 1) {
      modelBlocks[[i]] <- NULL
    }
  }






  #### Visualize ##############################################################


}

