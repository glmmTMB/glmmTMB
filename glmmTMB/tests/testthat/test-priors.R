stopifnot(require("testthat"),
          require("glmmTMB"))
## source("tests/testthat/setup_makeex.R")
cprior1 <- data.frame(prior = "normal(0, 3)",
                     class = "beta",
                     coef = "")

cprior2 <- data.frame(prior = "normal(0, 3)",
                     class = "beta",
                     coef = "period2")

cprior3 <- data.frame(prior = "normal(0, 3)",
                     class = "beta",
                     coef = 2)

cprior4 <- data.frame(prior = "gamma(1e8, 2)",
                     class = "ranef",
                     coef = "1|herd")

cprior5 <- data.frame(prior = "gamma(1e8, 2)",
                     class = "ranef_sd",
                     coef = "1|herd")

cprior6 <- data.frame(prior = "lkj(1)",
                     class = "ranef_cor",
                     coef = "1|herd")

gm0p1 <- update(gm0, priors = cprior1)
gm1p2 <- update(gm1, priors = cprior2)
gm1p3 <- update(gm1, priors = cprior3)

## try(update(gm1, priors = cprior4))

update(gm1, priors = cprior4)
update(gm1, priors = cprior5)

test_that("cor prior with no cor", {
    expect_warning(gm1p6 <- update(gm1, priors = cprior6))
    expect_equal(logLik(gm1p6), logLik(gm1))
})

get_prior_info <- function(fit) {
    pp <- fit$obj$env$data
    pvars <- grep("^prior", ls(pp), value = TRUE)
    mget(pvars, list2env(pp))
}

test_that("basic prior info", {
    expect_equal(get_prior_info(gm0p1),
                 list(prior_distrib = 0, prior_elend = 0, prior_elstart = 1,
                      prior_npar = 2,
                      prior_params = c(0, 3), prior_whichpar = 0))
})

test_that("prior printing", {
    cc <- capture.output(print(gm0p1))
    expect_equal(tail(cc[nzchar(cc)], 1), "Priors: fixef ~ normal(0, 3) ")
})

test_that("summary prior printing", {
    cc <- capture.output(print(summary(gm0p1)))
    expect_equal(trimws(tail(cc, 2)), c("Priors:", "fixef ~ normal(0, 3)"))
})

test_that("regularization example", {
    ## example from GH849
    set.seed(101)
    dd <- data.frame(d = ordered(rep(1:5, each = 10)))
    dd$resp <- ifelse(dd$d=="1", 0, rbinom(50, size = 1, prob = 0.5))
    m1 <- glmmTMB(resp ~ d, data = dd, family = binomial, prior = data.frame(prior = "normal(0,3)", class = "beta"))
    expect_equal(fixef(m1)$cond,
                 c(`(Intercept)` = -0.810093355244853, d.L = 1.42455679456297, 
                   d.Q = -2.48401389043977, d.C = 1.32524038308365, `d^4` = -0.231093396093939),
                 tolerance = 1e-5)
})

test_that("check for correct number of prior parameters", {
    priors3 <- data.frame(
        prior = "t(1,1)",
        class = "fixef"
    )
    expect_error(m3 <- glmmTMB(count ~ mined + (1 | site),
                  zi = ~mined,
                  family = poisson, data = Salamanders,
                  priors = priors3),
                 "was expecting")
})

test_that("print specific coefs for priors", {
    priors4 <- data.frame(
        prior = c("normal(0,2)", "t(1,1,10)"),
        class = c("fixef", "fixef"),
        coef = c("minedno", "DOP")
    )
    m4 <- glmmTMB(count ~ mined + DOP + (1 | site),
                  zi = ~mined,
                  family = poisson,
                  priors = priors4,
                  data = Salamanders
                  )
    cc <- capture.output(summary(m4))
    expect_true(any(grepl("fixef(minedno)", cc, fixed = TRUE)))
    expect_true(any(grepl("fixef(DOP)", cc, fixed = TRUE)))
})

