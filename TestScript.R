# Testing script

#### Define parameters for data pull ##########################################
NEONsites <- c("HOPB")
startdate <- "2018-01"
enddate <- "2018-12"

#### Pull raw data and K from NEON API (custom function) ######################
data <- request_NEON(NEONsites, startdate, enddate)
