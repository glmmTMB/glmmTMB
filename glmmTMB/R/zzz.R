## Startup code

# register emmeans methods dynamically
.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE))
        emmeans::.emm_register("glmmTMB", pkgname)
}
