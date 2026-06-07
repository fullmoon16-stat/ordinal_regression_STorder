# Data sources

## Mental health data

The mental health data are loaded directly from the R package `cgam` by using:

```r
library(cgam)
data(mental)
```
## Retinopathy data

The retinopathy data are loaded directly from the R package `catdata` by using:

```r
library(catdata)
data("retinopathy")
```

## Fine particulate matter (PM2.5) data

The PM2.5 data were constructed from air-quality records downloaded from public Seoul data sources. Since the original records form a time series, the final data set used for analysis was obtained by sampling observations to reduce temporal dependence.

The sampled data set used in the manuscript is provided as `pm25_seoul_sampled.csv`.
