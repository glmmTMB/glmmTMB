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
