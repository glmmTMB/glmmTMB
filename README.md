# glmmTMB

[![Build Status](https://travis-ci.org/glmmTMB/glmmTMB.svg?branch=master)](https://travis-ci.org/glmmTMB/glmmTMB)

`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It is currently pre-alpha or alpha software, but is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) and zero-inflation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

## Installation 

### From Github (source)

You can install `glmmTMB` via
```
devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
```
(this string denotes "Github user `glmmTMB`, repository `glmmTMB`, subdirectory `glmmTMB`"). If the install fails at the vignette-building step, try specifying `build_vignettes=FALSE` within the `install_github` call. Alternatively you can use `install_github()` from the `ghit` package, which has fewer dependencies. You'll need to have development tools (compilers etc.) installed first: `devtools::setup_rtools()` might help with that. Installing the source version will ensure that you get the very latest version of the package: since the package is in rapid development, that's a good idea. 

### From GH repository

If you have any problems with the source install (e.g. you're using MacOS or Windows and don't have the development tools), the following alternative incantation should work:
```
grepo <- "http://glmmtmb.github.io/glmmTMB/repos/"
options(repos=c(grepo,getOption("repos")))
install.packages("glmmTMB")
```
At present this repository is only updated manually, so if you find the package to be unavailable or out of date by this route, please contact the maintainers.


