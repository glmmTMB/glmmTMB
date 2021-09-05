## Retrieve from cache, or re-build, fitted objects that are used many times
## (rebuilt on every new installation *but* shared between e.g. different archs)
not_cran <- identical(Sys.getenv("NOT_CRAN"), "true")
up2date <- function(oldfit) {
  if (TMB:::isNullPointer(oldfit$obj$env$ADFun$ptr)) {
    obj <- oldfit$obj
    oldfit$obj <- TMB::MakeADFun(obj$env$data,
                                 obj$env$parameters, 
                                 map = obj$env$map,
                                 random = obj$env$random, 
                                 silent = obj$env$silent,
                                 DLL = "glmmTMB")
    oldfit$obj$env$last.par.best <- obj$env$last.par.best
  }
  oldfit
}
file_ok <- gt_load("test_data/models.rda")
if (!file_ok) {
    make_ex <- system.file("test_data", "make_ex.R", package="glmmTMB", mustWork = TRUE)
    cat("running fits to build test examples ...\n")
    source(make_ex, chdir = TRUE, echo = FALSE)
}
