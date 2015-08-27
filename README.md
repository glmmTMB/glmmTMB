# glmmTMB

`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It is currently pre-alpha or alpha software, but is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) and zero-inflation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

## Installation 

To install `glmmTMB`, first install TMB following [its installation instructions](https://github.com/kaskr/adcomp/blob/master/README.md), then install `glmmTMB` via
```
devtools::install_github("glmmTMB/glmmTMB",subdir="glmmTMB")
```
(you may need to install the `devtools` package, compiler tools, and other dependencies first).