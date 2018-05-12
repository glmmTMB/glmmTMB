# glmmTMB

[![Build Status](https://travis-ci.org/glmmTMB/glmmTMB.svg?branch=master)](https://travis-ci.org/glmmTMB/glmmTMB)
[![cran version](http://www.r-pkg.org/badges/version/glmmTMB)](https://cran.r-project.org/package=glmmTMB)
[![downloads](http://cranlogs.r-pkg.org/badges/glmmTMB)](http://cranlogs.r-pkg.org/badges/glmmTMB)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)


`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on [CppAD](https://www.coin-or.org/CppAD/) and [Eigen](eigen.tuxfamily.org/). It is relatively new, but is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) as well as model extensions such as zero-inflation, heteroscedasticity, and autocorrelation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

## Installation 

### From CRAN

`glmmTMB` is [on CRAN](https://CRAN.R-project.org/package=glmmTMB), so you can install the latest release version in the usual way, i.e. `install.packages("glmmTMB")` .

### Handling TMB/Matrix/glmmTMB mismatches

When loading `glmmTMB` you may encounter this message:

> Package version inconsistency detected.<br>
> TMB was built with Matrix version <xxxx>
> Current Matrix version is <yyyy>
> Please re-install 'TMB' from source or restore original 'Matrix' package

This occurs because you've updated the `Matrix` package to a newer version. Installing a new *binary* version of `TMB` probably won't help, because the binary package on CRAN have been built with the older version.

To re-install `TMB`, or to restore an older version of `Matrix`, you will need to have developer tools (compilers etc.) installed; these are not R packages, but additional packages and libraries for your operating system. You can try `devtools::dr_devtools()` to see if you have them already; if not, see the [RStudio devtools docs](https://www.rstudio.com/products/rpackages/devtools/) for links to download and install them.

At this point, if you're lucky, `install.packages("TMB",type="source")` will take care of everything.

If you're unlucky (e.g. you're using MacOS and originally installed R from a binary package), you may have some more work to do.

1. if you get errors about `library not found for -lgfortran` or `library not found for -lquadmath` you need to follow [these instructions](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/) to update your Fortran compilers
2. if you get errors about `unsupported option '-fopenmp'`you need to turn off OpenMP compilation by adding the line `SHLIB_OPENMP_CFLAGS=` to your `~/.R/Makevars` file. If you've already done step #1 this file will already exist; use a text editor to add the line above. Otherwise, you need to create it.

After updating your compilers (if necessary) and turning off OpenMP compilation, re-try the installation from source.

If you opt to restore an older version of `Matrix`, try `devtools::install_version("Matrix","<xxxx>")` (where `<xxxx>` is the version of `Matrix` shown in the original error message). This will also require installation from source.

If all else fails you can ask a maintainer to provide a binary version of the `TMB` package that works for your OS and `Matrix` version.

### From Github (source)

You can install the most recent development version of `glmmTMB` via
```
devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
```
(this string denotes "Github user `glmmTMB`, repository `glmmTMB`, subdirectory `glmmTMB`"). If the install fails at the vignette-building step, try specifying `build_vignettes=FALSE` within the `install_github` call. You'll need to have development tools (compilers etc.) installed: `devtools::dr_devtools()` and the [RStudio devtools docs](https://www.rstudio.com/products/rpackages/devtools/) should help. Installing the source version will ensure that you get the very latest version of the package.

## Where to ask questions

- `r-sig-mixed-models@r-project.org` for general questions about `glmmTMB` usage and mixed models
- https://github.com/glmmTMB/glmmTMB/issues for bug, infelicity, and wishlist reporting
- https://groups.google.com/forum/#!forum/tmb-users for TMB-specific questions
- maintainer e-mail only for urgent/private communications
