#' Estimate K at temperature, given temperature and
#' normalized K600 values. Internal function. See Wanninkhof 1992 and
#' Hall et al 2016 supplemental material
Kcor <- function(temp, K600) {
  K600 / (600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
