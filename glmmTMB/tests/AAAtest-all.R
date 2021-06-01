## Cache fitted objects that are used many times
## (rebuilt on every new installation *but* shared between e.g. different archs)
if (!file.exists(system.file("test_data", "models.rda", package="glmmTMB"))) {
  make_ex <- system.file("test_data", "make_ex.R", package="glmmTMB", mustWork=TRUE)
  source(make_ex, chdir=TRUE, echo=TRUE)
}
if(require("testthat")) {
    pkg <- "glmmTMB"
    require(pkg, character.only=TRUE)
    print(sessionInfo())
    test_check(pkg, reporter="summary")
    print(warnings()) # TODO? catch most of these by expect_warning(..)
} else {
    warnings("Package 'testthat' not available, cannot run unit tests for package",
	     sQuote(pkg))
}
