stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4")  ## for simulate
          )

context("glmmTMBControl")

## Some selected L1-distances between two fits
distFits <- function(fit1, fit2) {
    s1 <- summary(fit1)
    s2 <- summary(fit2)
    list(
        est_diff=max(abs((coef(s1)$cond - coef(s2)$cond)[,"Estimate"])),
        sd_diff=max(abs((coef(s1)$cond - coef(s2)$cond)[,"Std. Error"])),
        ll_diff=abs(logLik(fit1) - logLik(fit2))
    )
}

test_that("profile method", {

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

  myfit <- function(...) {
      glmmTMB(count ~ mined * spp + (1|site),
              zi = ~ mined,
              family = poisson,
              data = Salamanders,
              verbose = FALSE,
              control = glmmTMBControl(...))
  }

  # Record time and model
  capture_time_model <- function(..., FUN=myfit) {
    s1 <- system.time({
        model <- FUN(...)
    })
    return(list(model = model,
               time = s1))
  }


  m1 <- capture_time_model( parallel = 1 )
  dc <- parallel::detectCores()
  dc <- min(5,dc,getOption("mc.cores"))
  m2 <- capture_time_model( parallel = dc  )

  expect_true( all( distFits(m1[[1]], m2[[1]]) < c(1e-4, 1e-2, 1e-4) ) )

  if (FALSE) {
      if (dc >=2) {
          expect_true( m1$time[["elapsed"]] >= m2$time[["elapsed"]])
      }
  }
})

if (FALSE) {
    ## trying to do something bigger so the diffs are more obvious with
    ## only 2 cores ... ?
    set.seed(101)
    n <- 1e5
    nf <- 1e3
    dd <- data.frame(x=rnorm(n),f=factor(sample(nf,size=n,replace=TRUE)))
    dd$y <- simulate(~x+(x|f), data=dd,
                     family=gaussian,
                     newdata=dd,
                     newparams=list(beta=c("(Intercept)"=1,
                                           x=2),
                                    theta=c(1,1,2),
                                    sigma=1))[[1]]
    myfit2 <- function(...) {
        glmmTMB(y ~ x + (x|f),
                family = gaussian,
                data = dd,
                verbose = FALSE,
                control = glmmTMBControl(...))
    }
      
    capture_time_model(parallel=1,FUN=myfit2)
    capture_time_model(parallel=dc,FUN=myfit2)
}
