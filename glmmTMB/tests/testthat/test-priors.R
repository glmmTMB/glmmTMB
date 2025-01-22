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

test_that("print.priors_glmmTMB works without coef column #1014", {
  gdat <- readRDS(system.file("vignette_data", "gophertortoise.rds", package = "glmmTMB"))
  gmod_glmmTMB <- glmmTMB(
    shells ~ prev + offset(log(Area)) + factor(year) + (1 | Site),
    family = poisson,
    data = gdat
  )
  gprior <- data.frame(
    prior = "gamma(1e8, 2.5)",
    class = "ranef"
  )
  gmod_glmmTMB_p <- update(gmod_glmmTMB, priors = gprior)
  x <- summary(gmod_glmmTMB_p)
  out <- capture.output(print(x$priors))
  expect_identical(trimws(out), "ranef ~ gamma(1e+08, 2.5)")
  expect_no_error(capture.output(print(x)))
})

test_that("prior specs", {
    skip_on_cran()
    m2 <- glmmTMB(count ~ spp + mined + (1 | site),
                  zi = ~ spp + mined,
                  family = nbinom2, data = Salamanders
                  )
    prior6 <- data.frame(
        prior = c("normal(250,3)", "t(0,3,3)", "t(0,3,3)", "gamma(10,1)"),
        class = c("fixef", "fixef", "fixef_zi", "ranef_sd"),
        coef = c(1, 2, 1, 1)
    )
    g6p <- suppressWarnings(update(m2, prior = prior6)) ## NPD
    expect_equal(fixef(g6p)$cond[[1]], 248.6945, tolerance = 1e-5)
    
    prior7 <- data.frame(
        prior = c("normal(250,3)", "t(0,3,3)", "t(0,3,3)", "gamma(10,1)"),
        class = c("fixef", "fixef", "fixef_zi", "ranef_sd"),
        coef = c("(Intercept)", "minedno", "minedno", "(Intercept)")
    )
    expect_error(update(m2, prior = prior7), "can't match")

    prior8 <- prior7
    prior8$coef[4] <- "site"
    suppressWarnings(g8p <- update(m2, prior = prior8))

    prior8$coef[4] <- "1|site"
    suppressWarnings(g9p <- update(m2, prior = prior8))

    expect_equal(getME(g8p, "theta"), getME(g9p, "theta"))
})

## from GH#1146

rankdeficientdata <- data.frame(
  group = rep(c('g1', 'g2', 'g3'), c(10, 20, 10)),
  year = rep(c('y1', 'y2'), each = 20)
)
rankdeficientdata$y <- simulate_new(~group*year,
                  seed = 101,
                  newdata = rankdeficientdata,
                  newparams = list(beta=c(10,2,4,1), betadisp = log(3)))[[1]]

test_that("dropping priors for rank-def X matrix", {
    glmmTMB(y ~ group * year, data = rankdeficientdata,
            priors = data.frame(prior = c('normal(0, 10)'), class = c('fixef')))
})
