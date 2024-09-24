library(glmmTMB)
library(testthat)

data("sleepstudy", package = "lme4")
pp <- list(beta = c(280, 3),
           betadisp = 1,
           theta = c(3, 0.5, 2))

s1 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp)[[1]]

test_that("basic simulate_new",  {
    expect_equal(head(s1),
                 c(273.908009267299, 278.879666041338, 274.760596391018, 283.450108370691, 
                   286.469055736181, 290.145189139373))
})
                 

test_that("basic simulate_new, returning object", {
    s2 <- simulate_new( ~ Days + (Days|Subject),
                       seed = 101,
                       newdata = sleepstudy,
                       family = gaussian,
                       newparams = pp,
                       return_val = "object")
    expect_equal(fixef(s2)$cond,
                 c(`(Intercept)` = 280, Days = 3))
    expect_equal(c(VarCorr(s2)$cond$Subject),
                 rep(c(exp(3*2), exp(3)*exp(0.5)*2/sqrt(2^2+1), exp(0.5*2)), times = c(1, 2, 1)))
    ## SDs are exp(3), exp(0.5), cor is 2/(sqrt(2^2 +1)
    expect_equal(unname(head(getME(s2, "b"))),
                 c(-0.837259673343575,
                   6.24196658894249, 1.32388051566616, 
                   12.4287264928693, 0.829382894489882,
                   18.4190056),
                 tolerance = 1e-7)
})

s3 <- simulate_new( ~ Days + (Days|Subject),
                   newdata = sleepstudy,
                   seed = 101,
                   family = gaussian,
                   newparams = pp,
                   return_val = "pars")

test_that("simulate_new, returning pars", {
    ## FIXME: would be nice to have params in the canonical order
})

### now work on fixing particular elements ...
pp2 <- c(pp, list(b = rep(c(-1, 0, 1), length.out = 36)))

s4 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp2)[[1]]

test_that("simulate_new with specified b (values)", {
    ## should be different from random-b sim
    expect_false(isTRUE(all.equal(s1, s4)))
})

## need this for a downstream example
s5 <- simulate_new( ~ Days + (Days|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp2,
                   return_val = "pars")

test_that("simulate_new with specified b (pars)", {
    expect_identical(unname(s5[names(s5)=="b"]), pp2$b)
})                    

test_that("simulate_new with specified b (object)", {
    s6 <- simulate_new( ~ Days + (Days|Subject),
                       seed = 101,
                       newdata = sleepstudy,
                       family = gaussian,
                       newparams = pp2,
                       return_val = "object")
    expect_identical(unname(getME(s6, "b")),
                     ## correct, although names are wrong?
                     pp2$b)
})

pp3 <- c(pp, list(b = list("Days|Subject" = pp2$b)))

test_that("simulate_new with specified b (pars, b as list)", {
    s7 <- simulate_new( ~ Days + (Days|Subject),
                       seed = 101,
                       newdata = sleepstudy,
                       family = gaussian,
                       newparams = pp3,
                       return_val = "pars")
    expect_identical(unname(s5[names(s5)=="b"]), pp3$b[["Days|Subject"]])
})                    

pp4 <- list(beta = c(280),
           betadisp = 1,
           theta = c(-1, 1, 0))

s8 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp4)[[1]]

s9 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp4,
                   return_val = "pars")

test_that("simulate_new with two RE terms", {
    expect_equal(unname(head(s9)),
                 c(-0.119942121932298, 0.203239358640131, -0.248297964136576, 
                   0.078858438002708, 0.114325605998615, 0.431878061822679))

})

nb <- sum(names(s9) == "b")
pp5 <- c(pp4, list(b = rep(c(-1, 0, 1), length.out = nb)))

s10 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                    seed = 101,
                    newdata = sleepstudy,
                    family = gaussian,
                    newparams = pp5)[[1]]

test_that("simulate_new with two RE terms, fixed b", {
    s11 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                        seed = 101,
                        newdata = sleepstudy,
                        family = gaussian,
                        newparams = pp5, return_val = "pars")
    expect_equal(unname(s11[names(s11) == "b"]), pp5$b)
})

ns <- length(unique(sleepstudy$Subject))
pp6 <- c(pp4, list(b = list("0+factor(Days)|Subject" =
                                rep(c(-1, 0, 1), length.out = nb - ns))))


s12 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp6)[[1]]

test_that("simulate_new, b partially fixed", {
    expect_false(isTRUE(all.equal(s10, s12)))
})    

s13 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
                   seed = 101,
                   newdata = sleepstudy,
                   family = gaussian,
                   newparams = pp6,
                   return_val = "pars")

test_that("simulate_new, b partially fixed (pars)", {
    bvec <- unname(s13[names(s13) == "b"])
    expect_equal(head(bvec),
                 c(-0.119942121932298, 0.203239358640131, -0.248297964136576, 
                   0.078858438002708, 0.114325605998615, 0.431878061822679))
    expect_identical(bvec[-(1:ns)], pp6$b[[1]])
})

## s14 <- simulate_new( ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject),
##                     seed = 101,
##                    newdata = sleepstudy,
##                     family = gaussian,
##                    newparams = pp6,
##                     return_val = "object")

## Problem with these parameter entries:
## b 
## 3 
## Error in MakeADFun(data.tmb, parameters, map = mapArg, random = randomArg,  : 
##   Only numeric matrices, vectors and arrays can be interfaced
