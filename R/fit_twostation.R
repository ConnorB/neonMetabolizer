#' Model two-station stream metabolism
#'
#' @importFrom magrittr %>%
#'
#' @param data Dataframe containing formatted raw two station data, as
#'    returned by @request_NEON in $data
#'
#' @return
#'
#' @seealso
#'
#' @examples
#'
#' @export
fit_twostation <-function(data){
  #### Convert solarTime column to chron object ###############################
  # Mutate time column
  data <- data %>% mutate(date = lubridate::date(solarTime),
                          time = paste(lubridate::hour(solarTime),
                                       lubridate::minute(solarTime),
                                       lubridate::second(solarTime), sep = ":"))
  # Convert to chron object
  data$dtime <- chron(dates = as.character(data$date),
                      times = as.character(data$time),
                      format = c(dates = "y-m-d", times = "h:m:s"))
}

