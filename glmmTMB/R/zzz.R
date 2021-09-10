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
    checkDepPackageVersion(dep_pkg="TMB")
}
