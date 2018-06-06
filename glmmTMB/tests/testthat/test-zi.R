stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

## simulate something smaller than the Owls data set?


context("ZI models")
data(Owls)

test_that("zi", {
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

