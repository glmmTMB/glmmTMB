stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

## simulate something smaller than the Owls data set?


context("ZI models")
data(Owls)

test_that("zi", {
  skip_on_cran()
    ## Fit negative binomial model with "constant" Zero Inflation :
    owls_nb1 <<- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent +
                             (1|Nest)+offset(log(BroodSize)),
                         family = nbinom1(),
                         ziformula = ~1, data=Owls)
    owls_nb2 <<- update(owls_nb1,
                        ziformula = ~ FoodTreatment*SexParent + (1|Nest))
    owls_nb3 <<- update(owls_nb1,ziformula=~.)

    expect_equal(fixef(owls_nb2),
  structure(list(cond = structure(c(0.812028613585629, -0.342496105044418,
      -0.0751681324132088, 0.122484981295054),
      .Names = c("(Intercept)", "FoodTreatmentSatiated", "SexParentMale",
                 "FoodTreatmentSatiated:SexParentMale")),
      zi = structure(c(-2.20863281353936, 1.86779027553285, -0.825200772653965,
    0.451435813933449), .Names = c("(Intercept)", "FoodTreatmentSatiated",
     "SexParentMale", "FoodTreatmentSatiated:SexParentMale")),
    disp = structure(1.33089630005212, .Names = "(Intercept)")),
    .Names = c("cond", "zi", "disp"), class = "fixef.glmmTMB"),
  tolerance=1e-5)
    expect_equal(fixef(owls_nb2),fixef(owls_nb3))
})

test_that("zi beta and Gamma", {
    skip_on_cran()
    suppressWarnings(RNGversion("3.5.1"))
    set.seed(101)
    dd <- data.frame(yb=c(rbeta(100,shape1=2,shape2=1),rep(0,10)),
                     yg=c(rgamma(100,shape=1.5,rate=1),rep(0,10)))
    expect_error(glmmTMB(yb~1, data=dd, family=beta_family),
                 "y values must be")
    m1 <- glmmTMB(yb~1, data=dd, family=beta_family, zi=~1)
    expect_equal(unname(plogis(fixef(m1)[["zi"]])),1/11)
    expect_equal(unname(fixef(m1)[["cond"]]), 0.6211636, tolerance=1e-5)
    ## need *both* ziformula and family=ziGamma for gamma-hurdle
    expect_error(glmmTMB(yg~1, data=dd, family=Gamma),
                 "non-positive values not allowed")
    expect_error(glmmTMB(yg~1, zi=~1, data=dd, family=Gamma),
                 "non-positive values not allowed")
    expect_error(glmmTMB(yg~1, data=dd, family=ziGamma),
                 "non-positive values not allowed")
    m2 <- glmmTMB(yg~1, data=dd, family=ziGamma(link="log"), zi=~1)
    expect_equal(unname(plogis(fixef(m2)[["zi"]])),1/11)
    expect_equal(unname(fixef(m2)[["cond"]]), 0.3995267, tolerance=1e-5)
})
