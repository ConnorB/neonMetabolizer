#' Convert DateTime from UTC to local solar time
#' 
#' Convert DateTime from UTC to local solar time, which may be either apparent 
#' solar (perfect match between noon and solar zenith) or mean solar (exactly 24
#' hours between solar noons). This script is forked from the `streamMetabolizer` 
#' package (with some changes), and is reproduced here (instead of linking the 
#' package as a dependency) for forward compatability, as `streamMetabolizer` 
#' is no longer under active development, and dependencies `unitted` and 
#' `LakeMetabolizer` are not available for R versions >= 4.0.
#' 
#' @param date.time date-time values in POSIXct format and UTC timezone.
#' @param longitude numeric, in degrees, either positive and unitted ("degE" or 
#'   "degW") or with sign indicating direction (positive = East)
#' @param time.type character. "apparent solar", i.e. true solar time, is noon 
#'   when the sun is at its zenith. "mean solar" approximates apparent solar 
#'   time but with noons exactly 24 hours apart. Elsewhere in this package,
#'   variables named "solar.time" are mean solar time, whereas "app.solar.time"
#'   is apparent solar and "any.solar.time" is either.
#' @return a POSIXct object that says it's in tz="UTC" but that's actually in 
#'   solar time, with noon being very close to solar noon
#' @importFrom lubridate tz with_tz
#' @export
#' @references Yard, Bennett, Mietz, Coggins, Stevens, Hueftle, and Blinn. 2005.
#'   Influence of topographic complexity on solar insolation estimates for the 
#'   Colorado River, Grand Canyon, AZ. Ecological Modelling.
convert_UTC_to_solartime <- function(date.time, longitude, time.type = "mean solar"){
  # format checking - require tz=UTC and expected units
  if(class(date.time)[1] != "POSIXct")
    stop("expecting date.time as a POSIXct object")
  if(!(tz(date.time) %in% c("GMT","Etc/GMT-0","Etc/GMT+0","UTC"))) 
    stop("expecting tz=UTC")
  
  # Calculate mean.solar time, which approximates solar noon at clock noon to
  # within ~20 minutes
  longitude.UTC <- 0 # longitude of UTC time zone, units of degrees east
  adjustmentFactor <- 3.989
  time.adjustment <- adjustmentFactor * (longitude.UTC + abs(longitude)) # min per degree
  mean.solar <- date.time + as.difftime(time.adjustment, units = "mins")
  
  # Either return mean solar time or adjust to true (apparent) solar time
  if(time.type == "mean solar") {
    out <- mean.solar # Return mean solar time
  } else { # Return apparent solar time
    # Convert day of date and time to a decimal day of year, subtracting 1 to 
    # span dates between 0 and 364
    jday <- lubridate::day(mean.solar) + lubridate::hour(mean.solar)/24 + 
      lubridate::minute(mean.solar)/(60*24) - 1
    
    # Use the equation of time to compute the discrepancy between apparent and
    # mean solar time. E is in minutes.
    # Converting to radians with pi/180, equation of time as in Yard et al. 2005
    E <- 9.87 * sin(((2*360*(jday - 81))/365) * (pi/180)) -
      7.53 * cos(((360*(jday - 81))/365) * (pi/180)) -
      1.5 * sin(((360*(jday - 81))/365) * (pi/180)) # minutes
    out <- mean.solar + as.difftime(E, units = "mins")
  }
  return(out)
}