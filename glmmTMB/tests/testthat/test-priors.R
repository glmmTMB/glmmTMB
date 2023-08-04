stopifnot(require("testthat"),
          require("glmmTMB"))

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
                     class = "theta",
                     coef = 2)

gm0p1 <- update(gm0, priors = cprior1)
gm1p2 <- update(gm1, priors = cprior2)
gm1p3 <- update(gm1, priors = cprior3)

get_prior_info <- function(fit) {
    pp <- fit$obj$env$data
    pvars <- grep("^prior", ls(pp), value = TRUE)
    mget(pvars, list2env(pp))
}

test_that("basic prior info", {
    expect_equal(get_prior_info(gm0p1),
                 list(prior_distrib = 0, prior_elend = 0, prior_elstart = 0,
                      prior_params = c(0, 3), prior_whichpar = 0))
})

test_that("prior printing", {
    cc <- capture.output(print(gm0p1))
    expect_equal(tail(cc[nzchar(cc)], 1), "Priors: fixef ~ normal(0, 3) ")
})

test_that("summary prior printing", {
    cc <- capture.output(print(summary(gm0p1)))
    expect_equal(tail(cc, 2), c("Priors:", "fixef ~ normal(0, 3)"))
})

data("sleepstudy", package = "lme4")

fm1 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)
cprior <- data.frame(prior = "normal(0, 3)",
                     class = "beta",
                     coef = 2)

