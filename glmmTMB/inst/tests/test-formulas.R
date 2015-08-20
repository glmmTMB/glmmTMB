stopifnot(require("testthat"),
          require("glmmTMB"))

context("formula parsing")


nrt <- function(x) length(x$reTrmFormulas)

test_that("basic splitForm", {
  expect_equal(nrt(splitForm(y~(x+q))),0) ## reTrms part should be empty
  expect_equal(nrt(splitForm(y~(x+q)+(1|f))),1) 
  expect_equal(nrt(splitForm(y~(x+q)+us(1|f))),1) 
})