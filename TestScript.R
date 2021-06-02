# Testing script

#### Define parameters for data pull ##########################################
NEONsites <- c("HOPB")
startdate <- "2018-01"
enddate <- "2018-12"

#### Pull raw data and K from NEON API (custom function) ######################
data <- request_NEON(NEONsites, startdate, enddate)

# Save for testing
saveRDS(data, file = "data/Raw_HOPB_2018-01_2018-12.rds")
#data <- readRDS(file = "data/Raw_HOPB_2018-01_2018-12.rds")

#### Fit 2 station metabolism model to raw data ###############################
