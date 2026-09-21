stopifnot(require("testthat"),
          require("glmmTMB"))

## Some selected L1-distances between two fits
distFits <- function(fit1, fit2) {
    s1 <- summary(fit1)
    s2 <- summary(fit2)
    glmmTMB:::namedList(
        max(abs((coef(s1)$cond - coef(s2)$cond)[,"Estimate"])),
        max(abs((coef(s1)$cond - coef(s2)$cond)[,"Std. Error"])),
        abs(logLik(fit1) - logLik(fit2))
    )
}

test_that("profile method", {
  skip_on_cran()
    myfit <- function(...) {
        glmmTMB(count ~ mined * spp + (1|site),
                family = poisson,
                data = Salamanders,
                control = glmmTMBControl(...))
    }

    m1 <- myfit( profile=FALSE )
    m2 <- myfit( profile=TRUE  )

    expect_true( all( distFits(m1, m2) < c(1e-4, 1e-2, 1e-4) ) )

    ## ###########################################################

    myfit <- function(...) {
        glmmTMB(count ~ mined * spp + (1|site),
                zi = ~ (1 | spp),
                family = poisson,
                data = Salamanders,
                control = glmmTMBControl(...))
    }

    m1 <- myfit( profile=FALSE )
    m2 <- myfit( profile=TRUE  )

    expect_true( all( distFits(m1, m2) < c(1e-4, 1e-2, 1e-4) ) )
})

## data from GH #1317
set.seed(101)
d_1317 <- data.frame(y = rnorm(50))
d2_1317 <- data.frame(y = rnorm(50), x = rnorm(50))

test_that("profile=TRUE errors early with no free fixed effects (GH #1317)", {
    skip_on_cran()
    ctrl <- glmmTMBControl(profile = TRUE)
    fits <- list(
        no_fixef      = quote(glmmTMB(y ~ 0, data = d_1317, control = ctrl)),
        no_fixef_REML = quote(glmmTMB(y ~ 0, data = d_1317, control = ctrl,
                                      REML = TRUE)),
        all_mapped    = quote(glmmTMB(y ~ x, data = d2_1317, control = ctrl,
                                      start = list(beta = c(0, 0)),
                                      map = list(beta = factor(c(NA, NA))))),
        all_mapped_REML = quote(glmmTMB(y ~ x, data = d2_1317, control = ctrl,
                                        start = list(beta = c(0, 0)),
                                        map = list(beta = factor(c(NA, NA))),
                                        REML = TRUE)))
    for (nm in names(fits)) {
        expect_error(eval(fits[[nm]]), regexp = "profile", info = nm)
        expect_error(eval(fits[[nm]]), regexp = "free fixed-effect", info = nm)
    }
})

test_that("profile=TRUE works with a partially mapped beta (ML)", {
    skip_on_cran()
    expect_no_warning(
        m <- glmmTMB(y ~ x, data = d2_1317,
                     control = glmmTMBControl(profile = TRUE),
                     start = list(beta = c(0, 0)),
                     map = list(beta = factor(c(NA, 1))))
    )
    expect_equal(fixef(m)$cond[["(Intercept)"]], 0)
})

test_that("profile=TRUE guard ignores a fully mapped betazi", {
    skip_on_cran()
    ## conditional beta is free; only the zi coefficients are mapped
    ## ('$' partial matching of map$beta -> map$betazi would trip the guard)
    set.seed(102)
    d3 <- d2_1317
    d3$y <- rpois(50, exp(0.5 * d3$x))
    expect_no_error(
        glmmTMB(y ~ x, zi = ~ x, family = poisson, data = d3,
                control = glmmTMBControl(profile = TRUE),
                start = list(betazi = c(-3, 0)),
                map = list(betazi = factor(c(NA, NA))))
    )
})

test_that("profile=TRUE works with REML=TRUE", {
    skip_on_cran()
    cmp_reml <- function(label, ...) {
        m1 <- glmmTMB(..., REML = TRUE,
                      control = glmmTMBControl(profile = FALSE))
        m2 <- glmmTMB(..., REML = TRUE,
                      control = glmmTMBControl(profile = TRUE))
        expect_true( all( distFits(m1, m2) < c(1e-4, 1e-2, 1e-4) ),
                    info = label )
        expect_false( anyNA(vcov(m2, full = TRUE)), info = label )
    }
    cmp_reml("gaussian y ~ x", y ~ x, data = d2_1317)
    cmp_reml("salamanders", count ~ mined + (1|site),
             family = poisson, data = Salamanders)
    cmp_reml("salamanders zi", count ~ mined * spp + (1|site),
             zi = ~ (1|spp), family = poisson, data = Salamanders)

    ## Poisson with no random effects under REML: beta is integrated
    ## out and nothing else is estimated, so the rebuilt objective has
    ## an empty parameter vector and the Newton refinement must be skipped
    set.seed(103)
    d4 <- d2_1317
    d4$y <- rpois(50, exp(0.5 * d4$x))
    expect_no_warning(
        m2 <- glmmTMB(y ~ x, family = poisson, data = d4, REML = TRUE,
                      control = glmmTMBControl(profile = TRUE))
    )
    m1 <- glmmTMB(y ~ x, family = poisson, data = d4, REML = TRUE,
                  control = glmmTMBControl(profile = FALSE))
    expect_true( all( distFits(m1, m2) < c(1e-4, 1e-2, 1e-4) ),
                info = "poisson no-RE" )
})

test_that("whichNotRandom() drops the random-effect blocks, and beta when include_beta = TRUE", {
    nm <- c("beta", "b", "theta", "bzi", "betazi", "bdisp", "betadisp", "psi")
    expect_identical(whichNotRandom(nm), c(1L, 3L, 5L, 7L, 8L))
    expect_identical(whichNotRandom(nm, include_beta = TRUE),
                     c(3L, 5L, 7L, 8L))
    ## exact matching: "betazi"/"betadisp" are never dropped
    expect_identical(whichNotRandom(c("betazi", "betadisp"), include_beta = TRUE),
                     c(1L, 2L))
    expect_identical(whichNotRandom(character(0)), integer(0))
})

test_that("parallel regions", {

  skip_on_cran()

  myfit <- function(...) {
    glmmTMB(count ~ mined * spp + (1|site),
            family = poisson,
            data = Salamanders,
            verbose = FALSE,
            control = glmmTMBControl(...))
  }

  # Record time and model
  capture_time_model <- function(...) {
    start_time <- Sys.time()
    model <- myfit(...)
    end_time <- Sys.time()
    return(list(model = model,
                elapsed_time = end_time - start_time  ))
  }


  m1 <- capture_time_model( parallel = 1 )
  ## DON'T grab all cores - bad on large machines
  ## FIXME: check if parallel setting is persistent ???
  m2 <- capture_time_model( parallel = min(4, parallel::detectCores()  ))

  expect_true( all( distFits(m1[[1]], m2[[1]]) < c(1e-4, 1e-2, 1e-4) ) )

  # expect_true( m1[[2]] <= m2[[2]])
  

})

m1 <- glmmTMB(count ~ mined + (1|site), family = poisson, data = Salamanders)

test_that("autopar-only parallel", {
    m2 <- update(m1, control = glmmTMBControl(parallel = list(autopar = TRUE)))
    expect_equal(fixef(m1), fixef(m2))
})

get_parcore_output <- function(parallel_arg) {
    options(glmmTMB_openmp_debug = TRUE)
    cc <- capture.output(
        fit <- update(m1,
                      control = glmmTMBControl(parallel = parallel_arg))
    )
    options(glmmTMB_openmp_debug = FALSE)
    cc <- grep("NULL", cc, invert = TRUE, value = TRUE)
    n <- as.integer(gsub("\\D*(\\d+)\\D*", "\\1", cc[1], perl = TRUE))
    ap <- as.logical(gsub(".* autopar = (.*)", "\\1", cc[1], perl = TRUE))
    return(list(n = n, autopar = ap))
}

test_that("ncores/autopar argument handling", {
    ap0 <- getOption("glmmTMB.autopar", get_autopar())
    expect_equal(get_parcore_output(list(n = 2, autopar = TRUE)),
                 list(n=2L, autopar=TRUE))
    expect_equal(get_parcore_output(list(2, autopar = TRUE)),
                                    list(n=2L, autopar=TRUE))
    expect_equal(get_parcore_output(2),
                 list(n = 2L, autopar = ap0))
    expect_equal(get_parcore_output(list(n = 2)),
                 list(n = 2L, autopar = ap0))
})
