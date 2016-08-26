# glmmTMB

[![Build Status](https://travis-ci.org/glmmTMB/glmmTMB.svg?branch=master)](https://travis-ci.org/glmmTMB/glmmTMB)

`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It is currently pre-alpha or alpha software, but is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) and zero-inflation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

## Installation 

### From Github (source)

First install `TMB` from CRAN, then install `glmmTMB` via
```
devtools::install_github("glmmTMB/glmmTMB",subdir="glmmTMB")
```
You will need to make sure that the `devtools` package and development tools (compilers etc.) are installed first: `devtools::setup_rtools()` should help with this.

### From GH repository

If you are using MacOS or Windows and don't have the development tools, the following should work:
```
grepo <- "http://glmmtmb.github.io/glmmTMB/repos/"
options(repos=c(grepo,getOption("repos")))
install.packages("glmmTMB")
```
At present this repository is only updated manually, so if you find the package to be unavailable or out of date by this route, please contact the maintainers.


