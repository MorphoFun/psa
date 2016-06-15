# Phenotypic Selection Analysis tools (psa)

## Description
```psa``` is a package with functions to estimate the strength and mode and phenotypic selection. 

## Features
Some of the main features include:  
* Linear, quadratic, and correlational selection differentials, with option for bootstrapped standard errors
* Linear, quadratic, and correlational selection gradients, with option for bootstrapped standard errors
* Ordinary least-squares (OLS) regression method follows equations from Lande and Arnold (1983), so quadratic partial regression estimates and standard errors do NOT need to be doubled
* Data are standardardized to a mean of zero and unit variance, by default
* Options to compare selection coefficients and standard errors from different methods to evaluate robustness of estimates
* Option to evaluate how choice of traits influences estimates of selection gradients, when using the OLS regression

## Vignettes
Easy tutorials on using some of the ```psa``` functions are available through the following vignettes:
* How to estimate linear and nonlinear selection differentials, and comparing output from different methods
* How to estimate linear and nonliear selection gradients, and comparing output from different methods

## Installation
```psa``` is currently in a development version, and can be access via:

```
library(devtools)
install_github("MorphoFun/psa", build_vignettes = TRUE)
```

