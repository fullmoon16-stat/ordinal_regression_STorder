# Data sources

This file describes the data sets used in the empirical analyses in Section 5 of the manuscript. It provides the data sources, loading instructions, variable definitions, coding schemes, units, and preprocessing steps needed to understand the analysis data.


## Mental health data

The mental health data are loaded directly from the R package `cgam` by using:

```r
library(cgam)
data(mental)
```

The data set contains `n = 40` adult residents from Alachua County, Florida. The response variable is an ordered mental impairment status with four categories:

* `mental`: ordinal response variable, coded as 1 = well, 2 = mild symptom formation, 3 = moderate symptom formation, and 4 = impaired.

The covariates used in the manuscript are:

* `ses`: socioeconomic status, coded as 1 = high and 0 = low. This variable is denoted by SES in the manuscript.
* `life`: life-event index, an integer from 0 to 9 measuring the number and severity of important life events within the past three years. This variable is denoted by LEV in the manuscript.


## Retinopathy data

The retinopathy data are loaded directly from the R package `catdata` by using:

```r
library(catdata)
data("retinopathy")
```

The data set contains `n = 613` observations on persons with retinopathy-related information. The response variable is an ordered disease-severity variable:

* `RET`: ordinal response variable, coded as 1 = no retinopathy, 2 = nonproliferative retinopathy, and 3 = advanced/proliferative retinopathy or blindness.

The covariates used in the manuscript are:

* `DIAB`: diabetes duration in years.
* `GH`: glycosylated hemoglobin, measured in percent.
* `BP`: diastolic blood pressure, measured in mmHg.
* `SM`: smoking status, coded as 1 = smoker and 0 = non-smoker.

The continuous covariates `DIAB`, `GH`, and `BP` were standardized using the `scale()` function in R before model fitting.


## Fine particulate matter (PM2.5) data

The PM2.5 data set used in the manuscript was constructed from publicly available Seoul air-quality records and meteorological data. The final data set used for analysis is provided as `pm25_seoul_sampled.csv`.

Raw data sources:

* Seoul Open Data Plaza, daily average air-quality information by period:
  `https://data.seoul.go.kr/dataList/OA-2220/S/1/datasetView.do`

* Seoul Air Quality Information, daily average statistics by period:
  `https://cleanair.seoul.go.kr/statistics/dayAverage`

* Korea Meteorological Administration Open MET Data Portal, climate statistics by condition:
  `https://data.kma.go.kr/climate/RankState/selectRankStatisticsDivisionList.do`

The original air-quality records cover daily observations for the 25 districts of Seoul from January 2014 to May 1, 2025. District-level daily PM2.5 means were averaged to obtain one citywide daily PM2.5 value for each calendar day.

The response variable used in the analysis is:

* `grade`: ordered PM2.5 category, coded according to Korean air-quality forecast standards:

  * 1 = good, 0--15 $\mu\mathrm{g}/\mathrm{m}^3$,
  * 2 = moderate, $\mu\mathrm{g}/\mathrm{m}^3$,
  * 3 = bad, 36--75 $\mu\mathrm{g}/\mathrm{m}^3$,
  * 4 = very bad, above 75 $\mu\mathrm{g}/\mathrm{m}^3$.

The covariates used in the manuscript are:

* `Temp`: daily average temperature, measured in $^\circ\mathrm{C}$.
* `Humid`: humidity, measured in $\% rh$.
* `Wind`: wind speed, measured in $\mathrm{m}/\mathrm{s}$.
* `Press`: atmospheric pressure, measured in $\mathrm{hPa}$.
* `Rain`: rainfall indicator, coded as 0 = no rain and 1 = rainy.
* `Season`: seasonal indicator, coded as 0 = summer/fall and 1 = winter/spring.

Since the original records form a time series, the final data set used in the manuscript was obtained by random sampling with a minimum gap of 7 days between any two selected dates. The sampling was performed using seed 1534. Although the target sample size was set to 500, the minimum-gap constraint yielded 456 observations. The final sampled data set is provided as `pm25_seoul_sampled.csv`.

The continuous covariates were standardized before model fitting.
