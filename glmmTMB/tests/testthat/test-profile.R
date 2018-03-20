stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

pp <- NULL
fm2 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)

test_that("profiling", {
    pp <<- profile(fm2,npts=3, level_max=0.9)
    expect_equal(dim(pp),c(83,3))
    expect_equal(colnames(pp),c(".par",".focal","value"))
    ## FIXME: better naming ...
    expect_equal(levels(pp$.par),
                 c("(Intercept)", "Days", "d~(Intercept)", "theta_Days|Subject.1", 
                   "theta_Days|Subject.2", "theta_Days|Subject.3"))
})

test_that("confint",
          {
              expect_equal(confint(pp,level=0.9),
               structure(c(234.721869361326, 6.68846748689838, 6.22267801952836, 
                         2.54687304042687, 1.25468019450172, -0.67914339638892, 268.074294402097, 
                        14.2429229554084, 6.77182694221881, 3.71999315079224, 2.25548446321563, 
                        1.28181590832892), .Dim = c(6L, 2L), .Dimnames = list(c("(Intercept)", 
                        "Days", "d~(Intercept)", "theta_Days|Subject.1", "theta_Days|Subject.2", 
                        "theta_Days|Subject.3"), c("lwr", "upr"))),
               tolerance=1e-5)
          })
