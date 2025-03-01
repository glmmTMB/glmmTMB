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
