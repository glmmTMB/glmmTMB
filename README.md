# glmmTMB

[![Build Status](https://travis-ci.org/glmmTMB/glmmTMB.svg?branch=master)](https://travis-ci.org/glmmTMB/glmmTMB)
[![cran version](http://www.r-pkg.org/badges/version/glmmTMB)](https://cran.r-project.org/package=glmmTMB)
[![downloads](http://cranlogs.r-pkg.org/badges/glmmTMB)](http://cranlogs.r-pkg.org/badges/glmmTMB)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)
[![Research software impact](http://depsy.org/api/package/cran/glmmTMB/badge.svg)](http://depsy.org/package/r/glmmTMB)


`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It is currently pre-alpha or alpha software, but is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) and zero-inflation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

## Installation 

### From CRAN

`glmmTMB` is [on CRAN](https://CRAN.R-project.org/package=glmmTMB), so you can install the latest release version in the usual way, i.e. `install.packages("glmmTMB")`

### From Github (source)

You can install the very latest/development version of `glmmTMB` via
```
devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
```
(this string denotes "Github user `glmmTMB`, repository `glmmTMB`, subdirectory `glmmTMB`"). If the install fails at the vignette-building step, try specifying `build_vignettes=FALSE` within the `install_github` call. Alternatively you can use `install_github()` from the `ghit` package, which has fewer dependencies. You'll need to have development tools (compilers etc.) installed first: `devtools::setup_rtools()` might help with that. Installing the source version will ensure that you get the very latest version of the package: since the package is in rapid development, that's a good idea. 

