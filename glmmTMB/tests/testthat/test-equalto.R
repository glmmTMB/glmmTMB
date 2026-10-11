stopifnot(require("testthat"),
          require("glmmTMB"))

## ------ test with simulated data

test_that("equalto vs map-start on simulated data", {

  # simulate data for a multilevel meta-analysis 
  k.studies <- 10 #number of studies
  study <- rep(seq_len(k.studies), times=5) #study id - set 5 effect sizes per study
  k <- length(study) #total number of effect sizes
  id <- seq_len(k) #effect size id (across all studies) 
  es.id <- unlist(lapply(5, seq_len)) #id of within-study effect sizes
  b0 <- 0.2 #fixed effect coeff true value
  sigma2.u <- 0.2 #study level random effect
  sigma2.m <- 0.3 #obs level random effect
  set.seed(123); vi <- rbeta(k, 2, 20) # simulate sampling errors variances
  set.seed(123); u <- rnorm(k.studies, 0, sqrt(sigma2.u))[study] # simulate study level variance
  set.seed(123); m <- rnorm(k, 0, sqrt(sigma2.m)) # simulate obs level variance
  set.seed(123); e <- rnorm(k, 0, sqrt(vi))[study] #simulate samplign errors - without within-study correlation (just a diag VCV)
  y <- b0 + u + m + e #compute y
  
  # sim dataset
  dat <- data.frame(y = y, vi = vi, study = study, id = factor(id), es.id = es.id, g = 1)
  
  # fit with equalto
  V <- diag(dat$vi)
  fit1 <- glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, V),
                  data=dat,
                  REML=TRUE)
  # fit with map and start (fixing sampling errors with disp)
  fit2 <- glmmTMB(y ~ 1 + (1|study) + (1|id),
                     dispformula = ~ 0 + id,
                     map = list(betadisp = factor(rep(NA, nrow(dat)))), 
                     start = list(betadisp = log(sqrt(dat$vi))), 
                     data=dat,
                     REML=TRUE)
  
  expect_equal(c(VarCorr(fit1)$cond[[1]]),
               c(VarCorr(fit2)$cond[[1]]))
  
  expect_equal(c(sigma(fit1)^2),
               c(VarCorr(fit2)$cond[[2]][1]))
  
  # check these are equivalent 
  mod1 <- glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, diag(vi)), data=dat, REML=TRUE)
  mod2 <- glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, V), data=dat, REML=TRUE)
  expect_equal(mod1$fit$par, mod2$fit$par)
  expect_equal(mod1$sdr$gradient.fixed, mod1$sdr$gradient.fixed)

})

### checking error messages for input matrix 
# #-expect an error
# mod0 <- glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, vi), data=dat) 
# #-expect an error
# vix <- vi[1:45]
# modx <- glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, vix), data=dat)
# #-expect an error: one column matrix
# Vbad <- matrix(1:5)
# glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, Vbad), data=dat)
# #-expect an error: non-symmtric
# Vt <- matrix(1:5, 1:2)
# glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, Vt), data=dat)
# #-expect an error: input character vector
# a <- as.vector(c("a", "b", "c"))
# glmmTMB(y ~ 1 + (1|study) + equalto(0 + id|g, a), data=dat)


## ------ test comparing output with metafor with example dataset
test_that("compare glmmTMB equalto with metafor rma.mv", {
  
  skip_if_not_installed("metafor")
  suppressMessages(require(metafor))
  
  dat <- dat.assink2016
  dat$g <- 1
  dat$id <- as.factor(dat$id)
  dat$study <- as.factor(dat$study)
  
  # assume that the effect sizes within studies are correlated with rho=0.6
  V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)
  
  options(na.action = "na.omit") #metafor global options default
  fit.rma <- metafor::rma.mv(yi, V,
                    random = list(~ 1 | study, ~ 1 | id),
                    control=list(REMLf=FALSE), #to get the exact same likelihood as glmmTMB under REML
                    data=dat)
  
  fit.tmb <- glmmTMB(yi ~ 1 + (1 | study) + equalto(0 + id| g, V),
                     REML=TRUE,
                     data=dat)
  
  #overall mean est. 
  expect_equal(c(fit.rma$b[1]),
               c(summary(fit.tmb)$coefficients$cond[1]))
  
  #among study variance est.
  expect_equal(c(fit.rma$sigma2[1]),
               c(VarCorr(fit.tmb)$cond[[1]][1]),
               tolerance = 1e-6)

  #within study variance est.
  expect_equal(c(fit.rma$sigma2[2]),
               c(sigma(fit.tmb)^2),
               tolerance = 1e-6)
  
  #log-likelihood
  expect_equal(c(logLik(fit.rma)),
               c(logLik(fit.tmb)))
  
  #AIC
  expect_equal(c(AIC(fit.rma)),
               c(AIC(fit.tmb)))
})

