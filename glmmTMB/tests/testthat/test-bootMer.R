stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

context("bootMer")

test_that("bootMer works for Bernoulli responses", {
    fun <- function(x) predict(x)[1]
    Salamanders$pres <- as.numeric(Salamanders$count>0)
    m <- glmmTMB(pres ~ mined +(1|site), family=binomial, data=Salamanders)
    b <- lme4::bootMer(m, fun, nsim=2, seed=101)
    expect_true(var(c(b$t))>0)
})
