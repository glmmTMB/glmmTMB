## Startup code

## Work around RTMB (<= 2.0) registering an overly broad S4 solve() method
## that silently drops ... arguments: once RTMB's namespace is loaded, any
## solve(a, b, extra_arg) call *anywhere in the R session* breaks with
## "unused argument", because the method's signature ("num", "num.") matches
## any plain base numeric/matrix/array input, not just RTMB's own AD types.
## See https://github.com/kaskr/RTMB/issues/92
##
## RTMB's development version already fixes this (R/methods.R, commit
## "solve: add tol argument"): the method becomes
##   setMethod("solve", signature("num", "num."),
##             function(a, b, ...) base::solve(a, b, ...))
## Re-register the corrected method ourselves, here in .onLoad(): this runs
## after RTMB's namespace (imported via glmmTMB.R's @importFrom RTMB tag) is
## already loaded -- so RTMB's "num"/"num." classes exist -- but before
## glmmTMB's own namespace is locked, which is the only point at which
## setMethod(where=<default, our own namespace>) can succeed (namespaces are
## locked automatically once loading finishes, and by the time any other
## glmmTMB code runs, both glmmTMB's and RTMB's namespaces are already
## locked, so this cannot be done lazily on first use of the RTMB backend).
##
## Safe to remove once glmmTMB requires a fixed RTMB version.
fix_RTMB_solve <- function() {
    if (utils::packageVersion("RTMB") > "2.0") return(invisible(FALSE))
    methods::setMethod("solve", methods::signature("num", "num."),
                        function(a, b, ...) base::solve(a, b, ...))
    invisible(TRUE)
}

# register emmeans methods dynamically
.onLoad <- function(libname, pkgname) {
    fix_RTMB_solve()
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
