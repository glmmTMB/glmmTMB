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
    checkDepPackageVersion(dep_pkg="TMB")
}

.onUnload <- function(libpath) {
  library.dynam.unload("glmmTMB", libpath)
}
