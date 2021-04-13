# glmmTMB

[![cran version](http://www.r-pkg.org/badges/version/glmmTMB)](https://cran.r-project.org/package=glmmTMB)
[![downloads](http://cranlogs.r-pkg.org/badges/glmmTMB)](http://cranlogs.r-pkg.org/badges/glmmTMB)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)
[![R-CMD-check](https://github.com/glmmTMB/glmmTMB/workflows/R-CMD-check/badge.svg)](https://github.com/glmmTMB/glmmTMB/actions)

`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on [CppAD](https://www.coin-or.org/CppAD/) and [Eigen](eigen.tuxfamily.org/). It is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) as well as model extensions such as zero-inflation, heteroscedasticity, and autocorrelation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects models for the dispersion parameter.

## Where to ask questions

- the [r-sig-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for general questions about `glmmTMB` usage and mixed models (please subscribe to the list before posting)
- the [glmmTMB issues list](https://github.com/glmmTMB/glmmTMB/issues) for bug, infelicity, and wishlist reporting
- the [TMB users forum](https://groups.google.com/forum/#!forum/tmb-users) for TMB-specific questions
- maintainer e-mail only for urgent/private communications

Please do **not** cross-post, i.e. ask the same question in more than one forum, unless it's suggested that you have posted in the wrong place, or unless you receive total silence in one forum. In the latter case it may be better to send a reminder/"bump" message to the original forum; in any case you should mention in your new message where/when you've previously asked the question.

## Installation 

### From CRAN

`glmmTMB` is [on CRAN](https://CRAN.R-project.org/package=glmmTMB), so you can install the latest release version in the usual way, i.e. `install.packages("glmmTMB")` .

### From Github (source)

To get the very latest (development) version of `glmmTMB`, you can install it 
directly from the GitHub repository via
```
devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
```
(or using `remotes::install_github()` with the same argument); this string denotes "Github user `glmmTMB`, repository `glmmTMB`, subdirectory `glmmTMB`"). If the install fails at the vignette-building step, try specifying `build_vignettes=FALSE` within the `install_github` call. You will need to have development tools (compilers etc.) installed: `devtools::dr_devtools()` and the [RStudio devtools docs](https://www.rstudio.com/products/rpackages/devtools/) should help with that step. (Note of course that as well as being more up-to-date, the development version may contain new bugs or untested features.)

Recent versions of `glmmTMB` enable parallel (threaded) computations based on `OpenMP` (see the [parallel vignette](glmmTMB/vignettes/parallel.Rmd) for more information). OpenMP will be available automatically if your system supports it, but this may depend on the details of your operating system, compiler, compilation flags used when your R executable was built, etc.; in particular, see [here](https://github.com/Rdatatable/data.table/wiki/Installation#openmp-enabled-compiler-for-mac) for tips on enabling OpenMP for older (<= El Capitan/10.11.4) MacOS systems. (That page suggests using optimization level `-O3`, which [may cause problems for glmmTMB](https://github.com/glmmTMB/glmmTMB/issues/297).)

### Handling TMB/Matrix/glmmTMB mismatches

When loading `glmmTMB` you may encounter this message:

> Package version inconsistency detected.<br>
> TMB was built with Matrix version [xxxx]<br>
> Current Matrix version is [yyyy]<br>
> Please re-install 'TMB' from source or restore original 'Matrix' package

This may happen because you have installed a new version of `glmmTMB` from CRAN but haven't updated the `Matrix` package to its newest version. If this is the case, just use `update.packages()` (to update *all* of your packages) or `install.packages("Matrix")` (to install just the latest version of `Matrix` from CRAN).

Alternately (slightly more problematically), this may happen because you've updated the `Matrix` package to a newer version that was published on CRAN more recently than the latest CRAN version of `TMB`. Installing a new *binary* version of `TMB` from CRAN (i.e., via `update.packages()` or `install.packages("TMB")` on Windows or MacOS) probably won't help, because the binary package on CRAN will have been built with the older version.

**The easiest way to fix this problem** is to use the [checkpoint package]( https://CRAN.R-project.org/package=checkpoint) to revert your version of Matrix to the one that was available the last time the TMB package was updated on CRAN.

```
## load (installing if necessary) the checkpoint package
while (!require("checkpoint")) install.packages("checkpoint")
## retrieve build date of installed version of TMB
bd <- as.character(asDateBuilt(packageDescription("TMB",fields="Built")))
oldrepo <- getOption("repos")
setSnapshot(bd)
install.packages("Matrix")
options(repos=oldrepo) ## restore original repo
```

The only disadvantage to this approach is that your version of the `Matrix` package will be behind the version on CRAN; you might be missing out on some bug fixes or improvements, and eventually you may find that updates of other packages require the newest version of `Matrix`. (Also, if you accidentally update the `Matrix` package to the newest version, you'll have to redo this step.)

Alternatively, if you have development tools (compilers etc.) installed (see "installing from GitHub" above), `install.packages("TMB",type="source")` will take care of the problem, *if it works*.

If you're unlucky (e.g. you're using MacOS and originally installed R from a binary package), you may have some more work to do.

1. if you get errors about `library not found for -lgfortran` or `library not found for -lquadmath` you need to follow [these instructions](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/) to update your Fortran compilers.
2. if you get errors about `unsupported option '-fopenmp'`you need to turn off OpenMP compilation by adding the line `SHLIB_OPENMP_CFLAGS=` to your `~/.R/Makevars` file. If you've already done step #1 this file will already exist; use a text editor to add the line above. Otherwise, you need to create it.

After updating your compilers (if necessary) and turning off OpenMP compilation, re-try the installation from source.

If you opt to restore an older version of `Matrix`, try `devtools::install_version("Matrix","[xxxx]")`, where `[xxxx]` is the version of `Matrix` shown in the original error message. This will also require installation from source.

If all else fails you can ask a maintainer to provide a binary version of the `TMB` package that works for your OS and `Matrix` version.


