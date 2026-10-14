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


## Compare to asreml ----------------------------------------------------------------
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


  d <- agridat::kulkarni.sorghum.uniformity
  d$rowf <- factor(d$row); d$colf <- factor(d$col)
  d$yearf <- factor(d$year); d$one <- factor(1)

  fit_3sep <- glmmTMB(yield ~ yearf + sep(ar1(0 + rowf | one) %*% ar1(0 + colf | one) %*% ar1(0 + yearf | one) | one),
                        start = list(theta = rep(0, 6)),
                        map = list(theta = factor(c(1, 2, NA, 3, NA, 4))),
                        data = d, REML = TRUE)


# asreml3sep <- asreml(yield ~ yearf, random = ~ ar1(rowf):ar1(colf):ar1(yearf), data =d )
# asreml3sep <- update(asreml3sep)
# summary(asreml3sep)$varcomp  
  
  asreml_ref_3 <- list(rowcolyear_sigma2  = 737.1766373,   # rowf:colf:yearf 
                         rho_row = 0.9409499,   # rowf:colf:yearf!rowf!cor
                         rho_col = 0.2214698,   # rowf:colf:yearf!colf!cor
                         rho_year = 0.1800972  # rowf:colf:yearf!yearf!cor
                        )
  
  test_that("sep(ar1 %*% ar1) reproduces asreml estimates", {
    est <- VarCorr(fit_3sep)$cond
    expect_equal(unname(attr(est$one.sep.ar1.one., "stddev")[1]),  sqrt(asreml_ref_3$rowcolyear_sigma2),  tolerance = 1e-3)   
    expect_equal(unname(attr(est$one.sep.ar1.one., "correlation")[2, 1]), asreml_ref_3$rho_row, tolerance = 1e-3)     
    expect_equal(unname(attr(est$one.sep.ar1.one..1, "correlation")[2, 1]), asreml_ref_3$rho_col, tolerance = 1e-3)     
    expect_equal(unname(attr(est$one.sep.ar1.one..2, "correlation")[2, 1]), asreml_ref_3$rho_year, tolerance = 1e-3)   
  })
  
}
  
## ou and exp covstruc - times and dist components ----------------------------------------------------------------
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

## One test per covariance structure. 
## Each checks that sep(S %*% homdiag(group)) (homdiag log-sd is fixed at 0.)
## gives the same log-likelihood as the standard term S(... | group): an independent S per group. 
## Datasets (nlme): Orthodont (27 children x 4 ages) and Wheat2 (field trial).
## To do homtoep, toep

theta_est <- function(fit) {
  pf <- fit$sdr$par.fixed
  unname(pf[names(pf) == "theta"])
}

expect_sep_matches <- function(ref, term, lhs, data, grp, dispformula = ~1, start0 = 0) {
  n <- length(theta_est(ref))
  fit <- glmmTMB(
    as.formula(sprintf("%s + sep(%s %%*%% homdiag(0 + %s | one) | one)", lhs, term, grp)),
    dispformula = dispformula,
    start = list(theta = c(rep(start0, n), 0)),
    map   = list(theta = factor(c(seq_len(n), NA))),
    data  = data)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-4)
}

data(Orthodont, package = "nlme")
od <- as.data.frame(Orthodont)
od$Subject <- factor(od$Subject, ordered = FALSE)
od$agef  <- factor(od$age)          # 4 levels: factor structures
od$agenf <- numFactor(od$age)       # coordinates: ou / exp / gau
od$one   <- factor(1)

data(Wheat2, package = "nlme")
w <- as.data.frame(Wheat2)
w$latf <- factor(w$latitude)
w$lonf <- numFactor(w$longitude)
w$one  <- factor(1)

## Factor structures: dispformula = ~0 (one value per child and age, so no separate residual)
test_that("sep: diag", {
  ref <- glmmTMB(distance ~ age + diag(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "diag(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})
test_that("sep: homdiag", {
  ref <- glmmTMB(distance ~ age + homdiag(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "homdiag(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})
test_that("sep: homcs", {
  ref <- glmmTMB(distance ~ age + homcs(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "homcs(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})
test_that("sep: us", {
  ref <- glmmTMB(distance ~ age + us(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "us(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})
test_that("sep: cs", {
  ref <- glmmTMB(distance ~ age + cs(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "cs(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})

test_that("sep: ar1", {
  ref <- glmmTMB(distance ~ age + ar1(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "ar1(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})

test_that("sep: hetar1", {
  ref <- glmmTMB(distance ~ age + hetar1(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "hetar1(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})
test_that("sep: homcs", {
  ref <- glmmTMB(distance ~ age + homcs(0 + agef | Subject), dispformula = ~0, data = od)
  expect_sep_matches(ref, "homcs(0 + agef | one)", "distance ~ age", od, "Subject", ~0)
})
test_that("sep: rr", {
  ref <- glmmTMB(distance ~ age + rr(0 + agef | Subject, d = 2), data = od)
  expect_sep_matches(ref, "rr(0 + agef | one, d = 2)", "distance ~ age", od, "Subject",
                     start0 = 1)       # zero loadings are a degenerate start
})
## Coordinate structures 
test_that("sep: gau", {
  ref <- glmmTMB(distance ~ age + gau(0 + agenf | Subject), data = od)
  expect_sep_matches(ref, "gau(0 + agenf | one)", "distance ~ age", od, "Subject")
})
test_that("sep: mat", {      # Orthodont has too few ages to identify the shape; use Wheat2
  ref <- glmmTMB(yield ~ variety + mat(0 + lonf | latf), data = w)
  expect_sep_matches(ref, "mat(0 + lonf | one)", "yield ~ variety", w, "latf")
})
