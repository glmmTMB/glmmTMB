stopifnot(require("testthat"),
          require("glmmTMB"))

context("test edge cases")

test_that("profiling failure", {
    ## data from https://github.com/glmmTMB/glmmTMB/issues/399
    dd <- readRDS(system.file("test_data","IC_comp_data.rds", package="glmmTMB"))
    expect_warning(glmmTMB(
        ProbDiv ~ stdQlty + stdLaying + (1|Year) + (1|Site) + (1|PairID),
        family = "binomial",
        control=glmmTMBControl(profile = TRUE),
        data = dd),
        "a Newton step failed")
})


