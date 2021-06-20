stopifnot(require("testthat"),
          require("glmmTMB"))

## make sure tests don't run in parallel except where we want them to
op <- getOption("glmmTMB.cores", 1)
on.exit(options(glmmTMB.cores = op))
options(glmmTMB.cores = 1)

if (require(mvabund)) {

  data(spider)
  sppTot <- sort(colSums(spider$abund), decreasing = TRUE)
  tmp <- cbind(spider$abund, spider$x)
  tmp$id <- 1:nrow(tmp)
  spiderDat <- reshape(tmp,
                       idvar = "id",
                       timevar = "Species",
                       times =  colnames(spider$abund),
                       varying = list(colnames(spider$abund)),
                       v.names = "abund",
                       direction = "long")
  spiderDat_common <- subset(spiderDat, Species %in% names(sppTot)[1:4])


test_that("rr model fit", {
    ## Fit poison model with rr
    spider_p1 <<- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 1),
                         family = poisson,
                         data=spiderDat_common)
    spider_p2 <<- update(spider_p1,
                        control = glmmTMBControl(start_method = list(method = "res", jitter.sd = 0.2)))
    expect_equal(as.numeric(logLik(spider_p1)), c(-736.0022),
                 tolerance=1e-4)
    expect_equal(fixef(spider_p1),fixef(spider_p2),
                 tolerance=1e-4)
})

    test_that("rr works with symbolic d", {
    d <- 1
    spider_p1d <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = d),
                         family = poisson,
                         data=spiderDat_common)
    expect_equal(fixef(spider_p1d), fixef(spider_p1))
    expect_equal(VarCorr(spider_p1d), VarCorr(spider_p1))

    })

    test_that("rr error about un-eval'able d", {
              if (exists("junk")) rm(junk)
              expect_error(glmmTMB(abund ~ Species + rr(Species + 0|id, d = junk),
                                   family = poisson,
                                   data=spiderDat_common),
                           "can't evaluate reduced-rank dimension")
    })

    test_that("rr error about non-numeric d", {
        junk <- "junk"
        expect_error(glmmTMB(abund ~ Species + rr(Species + 0|id, d = junk),
                             family = poisson,
                             data=spiderDat_common),
                     "non-numeric value for reduced-rank dimension")
    })

    test_that("rr warning about parallel eval", {
        options(glmmTMB.cores = 2)
        expect_warning(
            glmmTMB(abund ~ Species + rr(Species + 0|id, d = 2),
                             family = poisson,
                    data=spiderDat_common)
            ## have to protect () in the regex ...
            , "rr\\(\\) not compatible with parallel execution")
        options(glmmTMB.cores = 1)
    })
    
    set.seed(101)
    n <- 1000
    ngrp <- 10
    dd <- data.frame(x=rnorm(n),y=rnorm(n),z=rnorm(n),
                     g1=factor(sample(ngrp,size=n,replace=TRUE)),
                     g2=factor(sample(ngrp,size=n,replace=TRUE)))
    beta <- 1; names(beta) <- "(Intercept)"
    theta <- rep(1,2*10)
    dd$w <- suppressMessages(simulate(~1 + (x+y+z|g1) + (x+y+z|g2),
                                      newdata=dd,
                                      family=gaussian,
                                      newparams=list(beta = beta,
                                                     theta=theta,sigma=1))[[1]])

    test_that("rr eigenvalues", {
        m1 <- glmmTMB(w ~ 1 + rr(x+y+z|g1,2), data=dd)
        eigenvalues <- zapsmall(eigen(VarCorr(m1)$cond$g1)$values)
        expect_equal(eigenvalues[3:4], c(0, 0))
        m2 <- glmmTMB(w ~ 1 + rr(x+y+z|g1,3)  + (x+y+z|g2), data=dd)
        eigenvalues <- zapsmall(eigen(VarCorr(m2)$cond$g1)$values)
        expect_equal(eigenvalues[4], 0)
    })

    ## FIXME: test, remove if unnecessary
    options(glmmTMB.control = op) ## just in case on.exit() is inappropriate?
} ## if require(mvabund)
