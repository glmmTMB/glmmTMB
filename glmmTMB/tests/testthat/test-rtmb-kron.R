## Tests for the RTMB kron() (Kronecker product) covariance structure

context("RTMB kron() covariance structure")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-6

set.seed(6001)
kron_data <- expand.grid(
  member = factor(1:2),
  time = factor(1:5),
  dyad = factor(1:25)
)
S_member <- matrix(c(1, 0.5, 0.5, 1.5), 2)
S_time <- 0.6^abs(outer(1:5, 1:5, "-"))
L <- t(chol(kronecker(S_time, S_member))) ## member varies fastest in the rows
b <- as.vector(L %*% matrix(rnorm(10 * 25), 10))
kron_data$x <- rnorm(nrow(kron_data))
kron_data$y <- 1 + 0.5 * kron_data$x + b + rnorm(nrow(kron_data), sd = 0.5)
## no time correlation
kron_data$y2 <- as.vector(t(chol(S_member)) %*% matrix(rnorm(250), 2)) +
  rnorm(nrow(kron_data), sd = 0.5)

fit_rtmb <- function(formula, data = kron_data) {
  glmmTMB(formula, data = data, control = glmmTMBControl(use_rtmb = TRUE))
}

fit_kron <- fit_rtmb(y ~ x + kron(us(0 + member) %x% ar1(0 + time) | dyad))

test_that("kron() matches an equivalent ordinary structure", {
  f1 <- fit_rtmb(y2 ~ kron(us(0 + member) %x% homdiag(0 + time) | dyad))
  f2 <- fit_rtmb(y2 ~ us(0 + member | dyad:time))
  expect_equal(logLik(f1), logLik(f2), tolerance = tol_logLik)
  expect_equal(unname(getME(f1, "theta")), unname(getME(f2, "theta")),
               tolerance = 1e-4)
})

test_that("kron() Gaussian likelihood matches the dense marginal likelihood", {
  vc <- VarCorr(fit_kron)$cond
  expect_equal(names(vc), c("dyad", "dyad.1"))
  ## as in kronecker(), the last (time) margin varies fastest
  B <- kronecker(vc[[1]], vc[[2]])
  Z <- as.matrix(getME(fit_kron, "Z"))
  V <- Z %*% kronecker(diag(25), B) %*% t(Z) +
    sigma(fit_kron)^2 * diag(nrow(kron_data))
  r <- kron_data$y - getME(fit_kron, "X") %*% fixef(fit_kron)$cond
  R <- chol(V)
  nll <- sum(log(diag(R))) + sum(backsolve(R, r, transpose = TRUE)^2) / 2 +
    nrow(kron_data) * log(2 * pi) / 2
  expect_equal(-as.numeric(logLik(fit_kron)), nll, tolerance = tol_logLik)
})

test_that("three-margin kron() matches ar1()", {
  dd <- expand.grid(a = factor(1:2), time = factor(1:4), b = factor(1:3),
                    g = factor(1:10))
  S <- kronecker(diag(3), kronecker(0.5^abs(outer(1:4, 1:4, "-")), diag(2)))
  dd$y <- as.vector(t(chol(S)) %*% matrix(rnorm(24 * 10), 24)) +
    rnorm(nrow(dd), sd = 0.5)
  f1 <- fit_rtmb(y ~ kron(homdiag(0 + a) %x% ar1(0 + time) %x% homdiag(0 + b) | g), dd)
  f2 <- fit_rtmb(y ~ ar1(0 + time | g:a:b), dd)
  expect_equal(logLik(f1), logLik(f2), tolerance = tol_logLik)
})

test_that("kron() works with formula(), simulate() and predict(newdata)", {
  expect_identical(deparse1(formula(fit_kron)),
                   "y ~ x + kron(us(0 + member) %x% ar1(0 + time) | dyad)")
  s <- simulate(fit_kron, nsim = 2, seed = 1)
  expect_equal(dim(s), c(nrow(kron_data), 2L))
  expect_true(all(is.finite(unlist(s))))
  expect_equal(predict(fit_kron, newdata = kron_data[1:10, ]),
               predict(fit_kron)[1:10])

  obj <- fit_kron$obj
  set_simcodes(obj, "zero")
  expect_equal(obj$simulate()$b, rep(0, 250))
  set_simcodes(obj, "fix")
  expect_equal(obj$simulate()$b, obj$report()$b)
  set_simcodes(obj, "random")
})

test_that("kron() gives clear errors", {
  expect_error(
    glmmTMB(y ~ kron(us(0 + member) %x% ar1(0 + time) | dyad),
            data = kron_data, control = glmmTMBControl(use_rtmb = FALSE)),
    "needs the RTMB back-end")
  expect_error(fit_rtmb(y ~ kron(rr(0 + member) %x% ar1(0 + time) | dyad)),
               "cannot be rr")
  expect_error(fit_rtmb(y ~ kron(us(member) %x% ar1(0 + time) | dyad)),
               "must look like")
})
