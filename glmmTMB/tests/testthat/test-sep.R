stopifnot(require("testthat"),
          require("glmmTMB"))
## Tests for the RTMB sep covariance structure
context("RTMB separable covariance structure")

skip_if_not_installed("RTMB")

####  load data ####
data("Salamanders", package = "glmmTMB")

Salamanders$one    <- factor(1)                           # single-level factor: component with no replication
spp_lvls <- levels(Salamanders$spp)
Salamanders$sample <- as.factor(Salamanders$sample)

namedC <- function(nm) {                          # any positive-definite C 
  C <- 0.5^abs(outer(seq_along(nm), seq_along(nm), "-"))
  dimnames(C) <- list(nm, nm)
  C
}

data("spider_long", package = "glmmTMB")
spider_long$one <- factor(1)

local_useRTMB(TRUE)

## formula parsing ----------------------------------------------------------------
test_that("sep terms parse into K component bars plus a replicate bar", {
  sp <- c(names(glmmTMB:::.valid_covstruct), "s")
  s1 <- reformulas::splitForm(y ~ x + sep(us(0 + a | g1) %*% ar1(0 + b | g2) | grp), specials = sp)
  s2 <- reformulas::splitForm(y ~ x + sep(us(0 + a | g1) %*% ar1(0 + b | g2)), specials = sp)
  expect_equal(s1$reTrmClasses, "sep")
  expect_length(s1$reTrmFormulas[[1]], 3)                          # 2 components + replicate bar
  expect_equal(s1$reTrmFormulas[[1]][[3]], quote(1 | grp))
  expect_equal(s2$reTrmFormulas[[1]][[3]], quote(1 | 1))           # no outer bar: dummy replicate
  expect_equal(s1$reTrmSep[[1]]$reTrmClasses, c("us", "ar1") )      # components only
})

tol = 1e-4
## Salamanders: two grouping factors, no outer bar 
## sep(diag %*% propto) is block diagonal over coefficients, so it equals propto x*spp
test_that("sep(diag %*% propto) equals propto on the interaction grouping", {
  Cs <- namedC(spp_lvls)                                                 # for 0 + spp
  Cx <- namedC(colnames(model.matrix(~ 0 + spp:DOP, Salamanders)))          # for 0 + spp:DOP (see note 1)
  f_sep <- count ~ spp + sep(diag(0 + DOP | site) %*% propto(0 + spp | sample, Cs))
  f_add <- count ~ spp + propto(0 + spp:DOP | site:sample, Cx)
  sep_fixed_diag_propto <-  glmmTMB(f_sep, Salamanders, family = poisson, start = list(theta = c(0 , 0)),
                                        map = list(theta = factor(c(1, NA))) )  # theta = (diag log sd, propto log lambda)
  propto_sal <- glmmTMB(f_add, Salamanders, family = poisson)
  expect_equal(as.numeric(logLik(sep_fixed_diag_propto)), as.numeric(logLik(propto_sal)), tolerance = tol)
  expect_equal(unname(fixef(sep_fixed_diag_propto)$cond), unname(fixef(propto_sal)$cond), tolerance = 1e-3)
})

test_that("sep(ar1 %*% homdiag | site) equals ar1 within each site:spp", {
  f_sep <- count ~ spp + sep(ar1(0 + sample | one) %*% homdiag(0 + spp | one) | site)
  f_ar1 <- count ~ spp + ar1(0 + sample | site:spp)
  sep_fixed_ar1 <- glmmTMB(f_sep, Salamanders,family = poisson,
                      start = list(theta = c(rep(0, 3))), map = list(theta = factor(c(1,2, NA))))       # theta = (ar1 log sd, ar1 rho, homdiag log sd)
  ar1_sal <- glmmTMB(f_ar1, Salamanders, family = poisson)
  expect_equal(as.numeric(logLik(sep_fixed_ar1)), as.numeric(logLik(ar1_sal)), tolerance = tol)
  expect_equal(unname(fixef(sep_fixed_ar1)$cond), unname(fixef(ar1_sal)$cond), tolerance = 1e-3)
})


if (require("agridat", quietly = TRUE)) {
  data(gilmour.serpentine)
  agDat <- gilmour.serpentine
  agDat$rowf <- factor(agDat$row) 
  agDat$colf <- factor(10*(agDat$col-8))
  agDat <- agDat[order(agDat$rowf, agDat$colf), ] # Sort order needed by asreml

  # library(asreml)
  # m1 <- asreml(yield ~ gen, data=agDat, resid = ~ar1(rowf):ar1(colf)); summary(m1)$varcomp  
  asreml_ref <- list(sigma2  = 19689.58935,   # rowf:colf!R
                     rho_row = 0.90418574,   # rowf:colf!rowf!cor
                     rho_col = 0.42463460    # rowf:colf!colf!cor
                    )
  
  agDat$g <- 1
  m1_glmmTMB <- glmmTMB(yield ~ gen + sep( ar1(0 + rowf| g) %*% ar1(0 + colf | g)), 
                          dispformula = ~0,
                          data = agDat, REML = TRUE, start = list(theta = c(0,0,0,0)),
                          map = list(theta = factor(c(1, 2, NA, 3))),
                          control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))  

  test_that("sep(ar1 %*% ar1) reproduces asreml estimates", {
    est <- VarCorr(m1_glmmTMB)$cond
    expect_equal(unname(attr(est$X.sep.ar1.g., "stddev")[1]),  sqrt(asreml_ref$sigma2),  tolerance = 1e-3)   
    expect_equal(unname(attr(est$X.sep.ar1.g., "correlation")[2, 1]), asreml_ref$rho_row, tolerance = 1e-3)     
    expect_equal(unname(attr(est$X.sep.ar1.g..1, "correlation")[2, 1]), asreml_ref$rho_col, tolerance = 1e-3)   
  })

  test_that("model fit does not depend on the order of the two factors in sep()", {
  m1b_glmmTMB <- glmmTMB(yield ~ gen + sep( ar1(0 + rowf| g) %*% ar1(0 + colf | g)), 
                          dispformula = ~0,
                          data = agDat, REML = TRUE,
                          start = list(theta = c(0,0,0,0)), map = list(theta = factor(c(1, 2, NA, 3))),
                          control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))  

  expect_equal(as.numeric(logLik(m1_glmmTMB)), as.numeric(logLik(m1b_glmmTMB)), tolerance = 1e-6)
  })  
}
  
  
if (require("nlme", quietly = TRUE)) {
  
test_that("sep(ou %*% homdiag) fits on Oxboys", {
  data(Oxboys, package = "nlme")
  Oxboys <- as.data.frame(Oxboys)
  Oxboys$Subject = factor(Oxboys$Subject, ordered = FALSE)
  Oxboys$agenf = numFactor(Oxboys$age)
  Oxboys$one = factor(1)

  ref <- glmmTMB(height ~ age + ou(0 + agenf | Subject), data = Oxboys)   # standard path
  sep_ou <- glmmTMB(height ~ age + sep(ou(0 + agenf | one) %*% homdiag(0 + Subject | one) | one),
                          start = list(theta = c(0,0,0)), map = list(theta = factor(c(1, 2, NA))),
                  data = Oxboys)
  
  expect_equal(as.numeric(logLik(sep_ou)), as.numeric(logLik(ref)), tolerance = 1e-4)
})
  
test_that("sep(exp %*% homdiag) matches exp(0 + lonf | latf) on Wheat2", {
  data(Wheat2, package = "nlme")
  w <- as.data.frame(Wheat2)
  w$lonf <- numFactor(w$longitude)
  w$latf <- factor(w$latitude)
  w$one  <- factor(1)

  ref <- glmmTMB(yield ~ variety + exp(0 + lonf | latf), data = w)
  sep_exp <- glmmTMB(
    yield ~ variety + sep(exp(0 + lonf | one) %*% homdiag(0 + latf | one)),
    start = list(theta = c(0, 0, 0)), map = list(theta = factor(c(1, 2, NA))),
    data = w)
  expect_equal(as.numeric(logLik(sep_exp)), as.numeric(logLik(ref)), tolerance = 1e-4)
})

}