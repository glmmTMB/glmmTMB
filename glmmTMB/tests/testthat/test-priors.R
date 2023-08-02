stopifnot(require("testthat"),
          require("glmmTMB"))

cprior <- data.frame(prior = "normal(0, 3)",
                     class = "beta",
                     coef = "")

rp1 <- list(prior_distrib = 0, prior_whichpar = 0,
            prior_element = NA_integer_, 
            prior_params = c(0, 3))

test_that("prior setup", {
    expect_equal(proc_priors(cprior), rp1)
    
    cprior2 <- cprior
    cprior2[["class"]] <- "fixef"
    expect_equal(proc_priors(cprior2), rp1)
})

gm0p <- update(gm0, priors = cprior)

test_that("prior printing", {
    cc <- capture.output(print(gm0p))
    expect_equal(tail(cc[nzchar(cc)], 1), "Priors: fixef ~ normal(0, 3) ")
})

test_that("summary prior printing", {
    cc <- capture.output(print(summary(gm0p)))
    expect_equal(tail(cc, 2), c("Priors:", "fixef ~ normal(0, 3)"))
})

data("sleepstudy", package = "lme4")

fm1 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)
cprior <- data.frame(prior = "normal(0, 3)",
                     class = "beta",
                     coef = 2)
## element by name, e.g. "Days" not allowed yet
proc_priors(cprior)
