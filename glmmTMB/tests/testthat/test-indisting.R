stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, package = "lme4")

## Create dyadic structure from sleepstudy.
## Treat pairs of subjects as dyads (9 dyads, 2 members each, 10 days).
sleepstudy2 <- sleepstudy
sleepstudy2$dyad   <- factor(rep(1:9, each = 20))
sleepstudy2$P1     <- as.integer(as.integer(sleepstudy2$Subject) %% 2 == 1)
sleepstudy2$P2     <- 1L - sleepstudy2$P1
sleepstudy2$Days_c <- sleepstudy2$Days - mean(sleepstudy2$Days)

## Fit models once and reuse across tests
fm_k1 <- suppressWarnings(
  glmmTMB(Reaction ~ Days +
            indisting(0 + P1 + P2 | dyad),
          data = sleepstudy2,
          REML = TRUE)
)

fm_k1_homcs <- suppressWarnings(
  glmmTMB(Reaction ~ Days +
            homcs(0 + P1 + P2 | dyad),
          data = sleepstudy2,
          REML = TRUE)
)

fm_k2 <- suppressWarnings(
  glmmTMB(Reaction ~ Days_c +
            indisting(0 + P1 + P2 + (P1+P2):Days_c | dyad) +
            homcs(0 + P1 + P2 | dyad:Days),
          data        = sleepstudy2,
          dispformula = ~0,
          REML        = TRUE)
)

## Extract variance-covariance components once
vc_k2   <- VarCorr(fm_k2)$cond[[1]]
sds_k2  <- attr(vc_k2, "stddev")
cors_k2 <- attr(vc_k2, "correlation")

test_that("indisting k=1: log-likelihood matches homcs", {
  expect_equal(logLik(fm_k1), logLik(fm_k1_homcs), tolerance = 1e-4)
})

test_that("indisting k=2 (4x4): symmetry constraints satisfied", {
  expect_equal(unname(sds_k2[1]), unname(sds_k2[2]), tolerance = 1e-4)
  expect_equal(unname(sds_k2[3]), unname(sds_k2[4]), tolerance = 1e-4)
  expect_equal(unname(cors_k2[1, 3]), unname(cors_k2[2, 4]), tolerance = 1e-4)
  expect_equal(unname(cors_k2[1, 4]), unname(cors_k2[2, 3]), tolerance = 1e-4)
})

test_that("indisting k=2 (4x4): theta count equals k*(k+1) + 2 homcs = 8", {
  expect_equal(length(glmmTMB::getME(fm_k2, "theta")), 8L)
})
