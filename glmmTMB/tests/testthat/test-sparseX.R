stopifnot(require("testthat"),
          require("glmmTMB"))

load(system.file("test_data", "models.rda", package="glmmTMB",
                 mustWork=TRUE))

context("sparse X models")

test_that("basic fits", {
    fm2S <- update(fm2,sparseX=c(cond=TRUE))
    expect_equal(fixef(fm2),fixef(fm2S))
    expect_equal(VarCorr(fm2),VarCorr(fm2S))
    expect_equal(dim(getME(fm2,"X")),dim(getME(fm2S,"X")))
    expect_equal(predict(fm2),predict(fm2S))
    nd <- sleepstudy[1,]
    expect_equal(predict(fm2, newdata=nd),predict(fm2S, newdata=nd))
})

