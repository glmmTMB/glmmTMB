## Retrieve from cache, or re-build, fitted objects that are used many times
## (rebuilt on every new installation *but* shared between e.g. different archs)
not_cran <- identical(Sys.getenv("NOT_CRAN"), "true")
saved_ex <- system.file("test_data", "models.rda", package="glmmTMB")
if (!not_cran && file.exists(saved_ex)) {
  load(saved_ex)
} else {
  make_ex <- system.file("test_data", "make_ex.R", package="glmmTMB", mustWork = TRUE)
  cat("running fits to build test examples ...\n")
  source(make_ex, chdir = TRUE, echo = FALSE)
}
