
# SimplexRegression

<!-- badges: start -->
<!-- badges: end -->

The goal of SimplexRegression is to ...

## Installation

You can install the development version of SimplexRegression from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("dudajustino/SimplexRegression")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SimplexRegression)
## Fit simplex regression model
model <- simplexreg(y ~ x1 + x2, data = your_data, 
                    link.mu = "plogit1", 
                    link.sigma2 = "log")

# Model summary
summary(model)
```

## Code of Conduct

Please note that the SimplexRegression project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
