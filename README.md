# neonMetabolizer
neonMetabolizer is an in-development R package that is intended to pull [National Ecological Observatory Network (NEON)](http://www.neonscience.org/) data and model two-station stream metabolism. This package is __not__ a NEON or U.S. Geological Survey (USGS) product and is __not__ supported, maintained, or endorsed by these agencies.

## Dependencies
NEONmetabolizer relies on [@NEONScience](https://github.com/NEONScience)'s [`NEON-utilities`](https://github.com/NEONScience/NEON-utilities), [`localPressureDO`](https://github.com/NEONScience/NEON-water-quality/localPressureDO), and [`reaRate`](https://github.com/NEONScience/NEON-reaeration/reaRate) packages. NEON-metabolizer also uses some functions from [@USGS-R](https://github.com/USGS-R)'s [`streamMetabolizer`](https://github.com/USGS-R/streamMetabolizer) package. 

## What does it do?

1. `request_NEON()`  
    - Pulls and formats water quality, temperature, PAR, nitrate, discharge, water chemistry, aquatic plant counts, and barometric pressure data products during the time period of choice into a single time-series dataframe
    - Pulls SF6 data to calculate K600 reaeration rates and the relationship between measured K and measured stream discharge
2. `clean_NEON()` 
    - Takes products returned from `request_NEON()` and checks for obviously erronious sensor data (i.e. things like discharge < 0), equal 15-minute breaks throughout the time series, and the amount of missing data 
    - Any missing data is filled using an ARIMA model (check out the documentation for [`auto.arima()`](https://www.rdocumentation.org/packages/forecast/versions/8.16/topics/auto.arima) for more info)
    - K600 and travel time between sensor stations are extrapolated for each time point using the relationship between field measurements and discharge 
    - Column names and formatting are changed to jive with `streamMetabolizer`'s requirements for metabolism modeling
3. Model metabolism
    - Functions for this step are in active development - see notes below!

## Disclaimer
This package is __not__ a product of NEON or USGS, and is __not__ supported, maintained, or endorsed by these agencies.
