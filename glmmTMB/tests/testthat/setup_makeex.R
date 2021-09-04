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
saved_ex <- system.file("test_data", "models.rda", package="glmmTMB")
cat(sprintf("%s exists=%s, not_cran=%s\n",saved_ex,file.exists(saved_ex), not_cran))
if (file.exists(saved_ex)) {
    cat("loading and updating test examples\n")
    L <- load(saved_ex)
    for (m in L) {
        if (inherits(m, "glmmTMB")) {
            cat(m,"\n")
            assign(m, up2date(get(m)))
        }
    }
} else {
  make_ex <- system.file("test_data", "make_ex.R", package="glmmTMB", mustWork = TRUE)
  cat("running fits to build test examples ...\n")
  source(make_ex, chdir = TRUE, echo = FALSE)
  
}
