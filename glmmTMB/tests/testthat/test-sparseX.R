stopifnot(require("testthat"),
          require("glmmTMB"))

## load(system.file("test_data", "models.rda", package="glmmTMB", mustWork=TRUE))

context("sparse X models")

test_that("basic fits", {
    fm2S <- update(fm2,sparseX=c(cond=TRUE))
    ## loosened from 1e-6 to 2e-6 for Solaris ...
    expect_equal(fixef(fm2), fixef(fm2S), tolerance=2e-6)
    expect_equal(VarCorr(fm2), VarCorr(fm2S),tolerance=1e-3)
    expect_equal(dim(getME(fm2,"X")), dim(getME(fm2S,"X")))
    expect_equal(predict(fm2), predict(fm2S), tolerance=1e-3)
    nd <- sleepstudy[1,]
    expect_equal(predict(fm2, newdata=nd),predict(fm2S, newdata=nd),
                 tolerance=1e-3)
})

test_that("back-compatibility", {
    x <- up2date(readRDS(system.file("test_data","oldfit.rds",
                             package="glmmTMB")))
    expect_is(VarCorr(x),"VarCorr.glmmTMB")
    expect_is(predict(x),"numeric")
})
