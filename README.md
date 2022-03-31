# NEONmetabolizer
NEONmetabolizer is a repo for R functions that can be used to pull [NEON](http://www.neonscience.org/) data from the API and model two-station stream metabolism. This is _not_ an official NEON package and is _not_ supported or endorsed by NEON.

## Dependencies
NEONmetabolizer relies on [@NEONScience](https://github.com/NEONScience)'s [`NEON-utilities`](https://github.com/NEONScience/NEON-utilities), [`localPressureDO`](https://github.com/NEONScience/NEON-water-quality/localPressureDO), and [`reaRate`](https://github.com/NEONScience/NEON-reaeration/reaRate). NEON-metabolizer also takes advantage of some functions from [@USGS-R](https://github.com/USGS-R)'s [`streamMetabolizer`](https://github.com/USGS-R/streamMetabolizer) R package. Most of the two-station modeling is based on the [Hall et al. 2016 Supplemental](https://www.doi.org/10.1007/s10021-015-9918-1) code, as implemented in [Kelly et al. 2021](https://doi.org/10.1029/2021JG006469).

## Example workflow

```r
# Define parameters for data pull
NEONsites <- c("HOPB")
startdate <- "2018-01" #YYYY-MM
enddate <- "2018-12"

# Pull raw data and K from NEON API
data <- request_NEON(NEONsites, startdate, enddate)
```

## Disclaimer
Function development is in early stages and subject to change; please feel free to fork and mod for your own purposes. This repo is _not_ a product of the National Ecological Observatory Network (NEON) or the U.S. Geological Survey (USGS), and is _not_ supported, maintained, or endorsed by these agencies.
