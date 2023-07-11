# glmmTMB

[![Build Status](https://travis-ci.org/glmmTMB/glmmTMB.svg?branch=master)](https://travis-ci.org/glmmTMB/glmmTMB)

`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) and zero-inflation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

The source code, issues lists, etc. are [on the main GitHub repository](https://github.com/glmmTMB/glmmTMB).

## Installation 

### Development version, from source 

To install `glmmTMB`, first install `TMB` from CRAN, then install `glmmTMB` via
```
remotes::install_github("glmmTMB/glmmTMB",subdir="glmmTMB")
```
(you may need to install the `remotes` package, compiler tools, and other dependencies first).

### Development version, from repository

A binary version of `glmmTMB` *may* be available for your operating system/version of R: try

```r
install.packages("glmmTMB", repos="https://glmmTMB.github.io/glmmTMB/repos")
```

If (1) this doesn't work, (2) you really need the development version, (3) you can't install from source as above, please contact the maintainers.

## Updated versions of TMB

Windows and MacOS binaries of TMB built with newer versions of the `Matrix` package *may* be available here as well. Try:

```r
install.packages("TMB", repos="https://glmmTMB.github.io/glmmTMB/repos")
```
