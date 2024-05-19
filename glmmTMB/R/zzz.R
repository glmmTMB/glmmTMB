## Startup code

# register emmeans methods dynamically
.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE)) {
        if (utils::packageVersion("emmeans") < "1.4") {
            warning("please install a newer version of emmeans (> 1.4)")
            return(NULL)
        }
        emmeans::.emm_register("glmmTMB", pkgname)
    }
    ## https://stackoverflow.com/questions/49056642/how-to-make-variable-available-to-namespace-at-loading-time/
    if (getRversion() < "4.4.0") {
        assign("%||%", function (x, y)  { if (is.null(x)) y else x }, envir = topenv())
    }
    check_dep_version(dep_pkg="TMB")
}

## https://github.com/lme4/lme4/issues/768
## https://github.com/kaskr/adcomp/issues/387
get_abi_version <- function() {
    if (utils::packageVersion("Matrix") < "1.6-2") return(numeric_version("0"))
    Matrix::Matrix.Version()[["abi"]]
}

.TMB.build.version <- packageVersion("TMB")

.onUnload <- function(libpath) {
  library.dynam.unload("glmmTMB", libpath)
}
