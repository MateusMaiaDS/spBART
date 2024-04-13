
# rspBART

<!-- badges: start -->
<!-- badges: end -->

The goal of rspBART is to ...

## Installation

You can install the development version of rspBART from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MateusMaiaDS/spBART")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rspBART)
data("airquality")
airquality <- airquality[complete.cases(airquality),]
airquality$Ozone <- (airquality$Ozone)^(1/3)
x.train <- airquality[,!(colnames(airquality) %in% "Ozone")]
y.train <- airquality$Ozone
spBART <- rspBART(x_train = x.train,nIknots = 10,
                  y_train = y.train,n_tree = 3,
                  x_test = x.train)
```

