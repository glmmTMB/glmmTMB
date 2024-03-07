# glmmTMB

<!-- badges: start -->
[![cran version](http://www.r-pkg.org/badges/version/glmmTMB)](https://cran.r-project.org/package=glmmTMB)
[![downloads](http://cranlogs.r-pkg.org/badges/glmmTMB)](http://cranlogs.r-pkg.org/badges/glmmTMB)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)](http://cranlogs.r-pkg.org/badges/grand-total/glmmTMB)
[![R-CMD-check](https://github.com/glmmTMB/glmmTMB/workflows/R-CMD-check/badge.svg)](https://github.com/glmmTMB/glmmTMB/actions)
<!-- badges: end -->

`glmmTMB` is an R package for fitting generalized linear mixed models (GLMMs) and extensions, built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It is intended to handle a wide range of statistical distributions (Gaussian, Poisson, binomial, negative binomial, Beta ...) and zero-inflation. Fixed and random effects models can be specified for the conditional and zero-inflated components of the model, as well as fixed effects for the dispersion parameter.

## Where to ask questions

- the [r-sig-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for general questions about `glmmTMB` usage and mixed models (please [subscribe to the list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) before posting)
- the [glmmTMB issues list](https://github.com/glmmTMB/glmmTMB/issues) for bug, infelicity, and wishlist reporting
- the [TMB users forum](https://groups.google.com/forum/#!forum/tmb-users) for TMB-specific questions
- maintainer e-mail only for urgent/private communications

Please do **not** cross-post, i.e. ask the same question in more than one forum, unless it's suggested that you have posted in the wrong place, or unless you receive total silence in one forum. In the latter case it may be better to send a reminder/"bump" message to the original forum; in any case you should mention in your new message where/when you've previously asked the question.

See [here](https://github.com/glmmTMB/glmmTMB) for the development site.

## Installation 

### Simple cases

- **from CRAN (release version)**: `install.packages("glmmTMB")`. (On Windows and MacOS this will install binary packages, by default: see below for installation from source and why you might want to do that.)
- **from GitHub (development version, from source)**: use `install.packages()` to install the `TMB` and `remotes` packages from CRAN, then `remotes::install_github("glmmTMB/glmmTMB/glmmTMB")`. If the install fails at the vignette-building step, try specifying `build_vignettes=FALSE` within the `install_github` call. You will need to have development tools (compilers etc.) installed: see [here](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) or [here](https://teuder.github.io/rcpp4everyone_en/020_install.html). (As well as being more up-to-date, the development version may contain new bugs or untested features!)
- **from GitHub (development version, from binary)**: A binary release of the development version *may* be available for your operating system/version of R from [the glmmTMB repository here](./repos/index.html). If a sufficient version isn't available and you are having insurmountable problems installing the current version of the package from source yourself, please contact the maintainers.

### Complications

### glmmTMB/TMB/Matrix mismatches

`glmmTMB` should be run on a system with the same version of `TMB` installed that was originally used to build the package; `TMB` and `Matrix` have a similar dependency. If you update your system with binary versions of the `TMB` or `Matrix` packages from CRAN that are *newer than the version of glmmTMB on CRAN*, you'll get warnings or error messages. 

You can:

- re-install the entire `Matrix` > `TMB` > `glmmTMB` stack *from source* (this is slight overkill, you might not need to re-install the whole stack, but it doesn't hurt):
    - Make sure that you have **development tools** installed. See [here](https://mac.r-project.org/tools/) for MacOS and [here](https://cran.r-project.org/bin/windows/Rtools/) for Windows (if you're on Linux, you probably know what you're doing with this ...) You will need both C++ and Fortran compilers (`gfortran`).
    - start a *clean* R/RStudio session (make sure no packages other than the base/core packages are loaded)
    - `install.packages("Matrix")`: this should install the latest version of `Matrix`
    - `install.packages("TMB", type = "source")` - this installs the latest version of `TMB` in a way that is *binary-compatible* with the latest Matrix
    - `install.packages("glmmTMB", type = "source")` - this installs the latest version of `glmmTMB` in a way that is *binary-compatible* with the latest TMB
- hope that updated binary versions are available here for your OS and R version. You can check `http://glmmtmb.github.io/glmmTMB/repos/bin/[OS]/contrib/[R_version]/PACKAGES`, where [OS] is "macosx" or "windows", and [R_version] is the *major* version of R you're using (e.g. 4.1). See "install development version from GitHub, binary" above. (Windows and MacOS binaries of TMB built with newer versions of the `Matrix` package *may* be available here as well. Try: `install.packages("TMB", repos="https://glmmTMB.github.io/glmmTMB/repos")`.)
- Use the [groundhog package](https://groundhogr.com/) to install binary versions of `Matrix`, `TMB`, and `glmmTMB` from a date when they were consistent with each other:

(**note**: this recipe is not well tested ... if you try it and encounter problems, please [post an issue](https://github.com/glmmTMB/glmmTMB/issues))

```r
## load (installing if necessary) the groundhog package
while (!require("groundhog", quietly=TRUE)) install.packages("groundhog")
## retrieve build date of installed version of glmmTMB
bd1 <- as.character(asDateBuilt(packageDescription("glmmTMB",fields="Built")))
groundhog.library("TMB", bd1)
bd2 <- as.character(asDateBuilt(packageDescription("TMB",fields="Built")))
groundhog.library("Matrix", bd2)
```
The only disadvantage to this approach is that your versions of `TMB` and `Matrix` will be behind the version on CRAN; you might be missing out on some bug fixes or improvements, and eventually you may find that updates of other packages require newer versions of these packages. (If you accidentally update the packages from CRAN, you'll have to redo this step.)
- Install older versions of `TMB` and/or `Matrix` *from CRAN, from source* using `remotes::install_version("[pkg]","[xxxx]")`, where `[pkg]` is TMB or Matrix and `[xxxx]` is the older package version referred to in the first error message you received.


### OpenMP compilation

*This section and the next may become outdated as new versions of operating systems are released; please notify the maintainers if you run into trouble.*

`glmmTMB` enables parallel (threaded) computations based on `OpenMP` (see the [parallel vignette](glmmTMB/vignettes/parallel.Rmd) for more information). OpenMP will be available automatically if your system supports it, but this may depend on the details of your operating system, compiler, compilation flags used when your R executable was built, etc.; in particular, see [here](https://github.com/Rdatatable/data.table/wiki/Installation#openmp-enabled-compiler-for-mac) for tips on enabling OpenMP for older (<= El Capitan/10.11.4) MacOS systems. (That page suggests using optimization level `-O3`, which [may cause problems for glmmTMB](https://github.com/glmmTMB/glmmTMB/issues/297).)

The maximum number of threads used defaults to 48; to increase this value when installing from source, you can use

```r
withr::with_makevars(c(PKG_CPPFLAGS="-DCPPAD_MAX_NUM_THREADS=128"), {
  remotes::install_github("glmmTMB/glmmTMB/glmmTMB")
}, assignment="+="
)
```
(or the equivalent for `remotes::install_cran()`).

### MacOS compilation issues

1. if you get errors about `library not found for -lgfortran` or `library not found for -lquadmath` you may need to follow [these instructions](https://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks--lgfortran-and--lquadmath-error/) to update your Fortran compilers.
2. if you get errors about `unsupported option '-fopenmp'`you need to turn off OpenMP compilation by adding the line `SHLIB_OPENMP_CFLAGS=` to your `~/.R/Makevars` file. If you've already done step #1 this file will already exist; use a text editor to add the line above. Otherwise, you need to create it.

After updating your compilers (if necessary) and turning off OpenMP compilation, re-try the installation from source.





