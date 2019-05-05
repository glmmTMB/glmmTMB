stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

test_that("error messages for user-spec start", {
    expect_error(
        glmmTMB(Reaction~Days+(Days|Subject), sleepstudy,
                start=list(beta=c(2))),
        "parameter vector length mismatch.*length\\(beta\\)==1, should be 2")
    expect_error(glmmTMB(Reaction~Days+(Days|Subject), sleepstudy,
                         start=list(junk=5)),
                 "unrecognized vector 'junk'")
})
