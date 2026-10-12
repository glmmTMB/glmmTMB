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

## fit model with one redundant par held at 1 (theta = 0 when theta = log(sd))
fit_fix <- function(f, data, fix, ...) {
  th <- glmmTMB(f, data, doFit = FALSE, ...)$parameters$theta
  th[fix] <- 0
  m <- factor(seq_along(th))
  m[fix] <- NA
  glmmTMB(f, data, start = list(theta = th), map = list(theta = m), se = FALSE, ...)
}

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
  sep_fixed_diag_propto <- fit_fix(f_sep, Salamanders, fix = 2, family = poisson)         # theta = (diag log sd, propto log lambda)
  propto_sal <- glmmTMB(f_add, Salamanders, family = poisson)
  expect_equal(as.numeric(logLik(sep_fixed_diag_propto)), as.numeric(logLik(propto_sal)), tolerance = tol)
  expect_equal(unname(fixef(sep_fixed_diag_propto)$cond), unname(fixef(propto_sal)$cond), tolerance = 1e-3)
})


test_that("sep(ar1 %*% homdiag | site) equals ar1 within each site:spp", {
  f_sep <- count ~ spp + sep(ar1(0 + sample | one) %*% homdiag(0 + spp | one) | site)
  f_ar1 <- count ~ spp + ar1(0 + sample | site:spp)
  sep_fixed_ar1 <- fit_fix(f_sep, Salamanders, fix = 3, family = poisson)       # theta = (ar1 log sd, ar1 rho, homdiag log sd)
  ar1_sal <- glmmTMB(f_ar1, Salamanders, family = poisson)
  expect_equal(as.numeric(logLik(sep_fixed_ar1)), as.numeric(logLik(ar1_sal)), tolerance = tol)
  expect_equal(unname(fixef(sep_fixed_ar1)$cond), unname(fixef(ar1_sal)$cond), tolerance = 1e-3)
})


