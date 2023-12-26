## Retrieve from cache, or re-build, fitted objects that are used many times
## (rebuilt on every new installation *but* shared between e.g. different archs)
not_cran <- identical(Sys.getenv("NOT_CRAN"), "true")
## CAREFUL: update_gauss_disp = TRUE only for models stored *before*
##  var -> SD reparameterization
## (store API information in model object?)
file_ok <- gt_load("test_data/models.rda", update_gauss_disp = TRUE)
if (!file_ok) {
    make_ex <- system.file("test_data", "make_ex.R", package="glmmTMB", mustWork = TRUE)
    cat("running fits to build test examples ...\n")
    source(make_ex, chdir = TRUE, echo = FALSE)
}
