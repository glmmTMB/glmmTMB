stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

## cat("ON CRAN:", testthat:::on_cran(), "\n")
context("bootMer")

fun <- function(x) predict(x)[1]

test_that("Bernoulli responses", {
  skip_on_cran()
    Salamanders$pres <- as.numeric(Salamanders$count>0)
    m <- glmmTMB(pres ~ mined +(1|site), family=binomial, data=Salamanders)
    b <- lme4::bootMer(m, fun, nsim=2, seed=101)
    expect_true(var(c(b$t))>0)
    expect_equal(suppressWarnings(c(confint(b))),
                 c(-1.579923,-1.250725),tolerance=1e-5)
})

test_that("binomial responses", {
    skip_on_cran()
    m <- glmmTMB(count ~ mined + (1|site), family=poisson, data=Salamanders)
    ss1 <- simulate(m,nsim=2,seed=101)
    b <- bootMer(m, fun, nsim=2, seed=101)
    expect_true(var(c(b$t))>0)
    expect_equal(suppressWarnings(c(confint(b))),
                 c(-0.7261239,-0.6921794),
                 tolerance=1e-5)
})
