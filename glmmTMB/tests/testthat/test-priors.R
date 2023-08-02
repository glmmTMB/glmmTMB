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
print(summary(gm0p))
## glmmTMB:::print.glmmTMB_prior(gm0p$modelInfo$priors, compact = TRUE)
