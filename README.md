# glmmTMB

[![Build Status](https://travis-ci.org/glmmTMB/glmmTMB.svg?branch=master)](https://travis-ci.org/glmmTMB/glmmTMB)
[![cran version](http://www.r-pkg.org/badges/version/glmmTMB)](https://cran.r-project.org/package=glmmTMB)
[![downloads](http://cranlogs.r-pkg.org/badges/glmmTMB)](http://cranlogs.r-pkg.org/badges/glmmTMB)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)


`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on [CppAD](https://www.coin-or.org/CppAD/) and [Eigen](eigen.tuxfamily.org/). It is relatively new, but is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) as well as model extensions such as zero-inflation, heteroscedasticity, and autocorrelation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

## Installation 

### From CRAN

`glmmTMB` is [on CRAN](https://CRAN.R-project.org/package=glmmTMB), so you can install the latest release version in the usual way, i.e. `install.packages("glmmTMB")` .

### From Github (source)

You can install the very latest/development version of `glmmTMB` via
```
devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
```
(this string denotes "Github user `glmmTMB`, repository `glmmTMB`, subdirectory `glmmTMB`"). If the install fails at the vignette-building step, try specifying `build_vignettes=FALSE` within the `install_github` call. Alternatively you can use `install_github()` from the `ghit` package, which has fewer dependencies. You'll need to have development tools (compilers etc.) installed first: `devtools::setup_rtools()` might help with that. Installing the source version will ensure that you get the very latest version of the package: since the package is in rapid development, that's a good idea. 

## Where to ask questions

- `r-sig-mixed-models@r-project.org` for general questions about `glmmTMB` usage and mixed models
- https://github.com/glmmTMB/glmmTMB/issues for bug, infelicity, and wishlist reporting
- https://groups.google.com/forum/#!forum/tmb-users for TMB-specific questions
- maintainer e-mail only for urgent/private communications