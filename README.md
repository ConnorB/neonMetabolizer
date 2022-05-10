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

## Notes on Two-Station Modeling

At the moment (May 2022), scripts for two-station modeling are in active development. 

I began by writing up two-station Bayesian modeling functions based on the [Hall et al. 2016 Supplemental](https://www.doi.org/10.1007/s10021-015-9918-1) code, as implemented in [Kelly et al. 2021](https://doi.org/10.1029/2021JG006469). This was working... _fine_, but the accuracy of model predictions wasn't always great, predicted K didn't always jive with measured K from NEON's SF6 reaeration numbers (even with a relatively constrained prior), and, of course, computation time was long. I found that the Bayesian results were requiring more quality checking and fine-combing, which again is... _fine_, but not the _most_ ideal if the goal is to model long time series at many sites.

I'm now leaning towards using a two-station MLE model for most instances. For this structure, I'm feeding the MLE a fixed daily K input based on the relationship between reaeration and discharge derived from NEON's SF6 tracer measurements, which (so far) seems to be a pretty strong relationship at many sites (yay!). As an initial test, I ran single-station MLEs' with fixed daily K using `metab_mle()` from [@USGS-R](https://github.com/USGS-R)'s [`streamMetabolizer`](https://github.com/USGS-R/streamMetabolizer), which approximated measured O2 surprisingly well (and were, of course, nice and quick, too).

I'm currently working on building out a forked copy of streamMetabolizer, available at [`michelleckelly/streamMetabolizer`](https://github.com/michelleckelly/streamMetabolizer), to include two-station MLE modeling within `metab_mle()` and it's associated functions. If you're reading this during JASM 2022, my current sticking point is that I've got an error in how I've specified the differential two-station MLE equation (within [`create_calc_dDOdt()`](https://github.com/michelleckelly/streamMetabolizer/blob/main/R/create_calc_dDOdt.R)), so my goal for early summer 2022 is to comb through that and iron everything out.

## Example workflow

```r
# Define parameters for data pull
NEONsites <- c("CUPE") # Rio Cupeyes, PR
startdate <- "2019-01" #YYYY-MM
enddate <- "2019-12"
plotPath <- "/data/plots" # Where you would like to save travel time plots?

# Pull raw NEON data (water quality, temperature, PAR, nitrate, 
# discharge, water chemistry, aquatic plant counts, 
# barometric pressure) and reaeration from NEON API
data <- request_NEON(NEONsites, startdate, enddate, plotPath)

# Clean NEON data
cleanedData <- clean_NEON(data = data$data, 
                          k600_clean = data$k600_clean,
                          k600_fit = data$k600_fit)
```

## Disclaimer
Package development is in early stages and subject to change; please feel free to fork and mod for your own purposes. Again, this package is __not__ a product of NEON or USGS, and is __not__ supported, maintained, or endorsed by these agencies. Right now, it's just supported by me, is very beta, and likely has bugs. :bug::bug:
