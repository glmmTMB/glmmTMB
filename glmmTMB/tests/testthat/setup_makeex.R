## Cache fitted objects that are used many times
## (rebuilt on every new installation *but* shared between e.g. different archs)
if (!file.exists(system.file("test_data", "models.rda", package="glmmTMB"))) {
  make_ex <- system.file("test_data", "make_ex.R", package="glmmTMB", mustWork=TRUE)
  cat("running fits to build test examples ...\n")
  source(make_ex, chdir=TRUE, echo = FALSE)
} else {
   load(system.file("test_data","models.rda",package="glmmTMB",
                    mustWork=TRUE))
}
