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

The PM2.5 data set used in the manuscript was constructed from publicly available Seoul air-quality records and meteorological data. The final data set used for analysis is provided as `pm25_seoul_sampled.csv`.

Raw data sources:

* Seoul Open Data Plaza, daily average air-quality information by period:
   `https://data.seoul.go.kr/dataList/OA-2220/S/1/datasetView.do`

* Seoul Air Quality Information, daily average statistics by period:
   `https://cleanair.seoul.go.kr/statistics/dayAverage`

* Korea Meteorological Administration Open MET Data Portal, climate statistics by condition:
   `https://data.kma.go.kr/climate/RankState/selectRankStatisticsDivisionList.do`

The original records cover daily observations for the 25 districts of Seoul from January 2014 to May 1, 2025. Since the original records form a time series, the final data set used in the manuscript was obtained by sampling observations to reduce temporal dependence.

