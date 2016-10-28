stopifnot(require("testthat"),
          require("glmmTMB"))

test_that("addForm", {
    expect_equal(addForm(y~x,~1,~z),y~x+1+z)
    expect_warning(addForm(y~x,z~1),
                   "discarding LHS")
})
