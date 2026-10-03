stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, package = "lme4")

## Create dyadic structure from sleepstudy.
## Treat pairs of subjects as dyads (9 dyads, 2 members each, 10 days).
sleepstudy2        <- sleepstudy
sleepstudy2$dyad   <- factor(rep(1:9, each = 20))
sleepstudy2$member <- factor(as.integer(as.integer(sleepstudy2$Subject) %% 2) + 1L)
sleepstudy2$Days_c <- sleepstudy2$Days - mean(sleepstudy2$Days)
sleepstudy2$Days_f <- factor(sleepstudy2$Days)

## Fit models once and reuse across tests
fm_k1 <- suppressWarnings(
  glmmTMB(Reaction ~ Days +
            indisting(0 + member | dyad),
          data        = sleepstudy2,
          dispformula = ~0,
          REML        = TRUE)
)

fm_k1_homcs <- suppressWarnings(
  glmmTMB(Reaction ~ Days +
            homcs(0 + member | dyad),
          data        = sleepstudy2,
          dispformula = ~0,
          REML        = TRUE)
)

fm_k2 <- suppressWarnings(
  glmmTMB(Reaction ~ Days_c +
            indisting(0 + member + member:Days_c | dyad) +
            homcs(0 + member | dyad:Days_f),
          data        = sleepstudy2,
          dispformula = ~0,
          REML        = TRUE)
)

## Extract variance-covariance components once
vc_k2   <- VarCorr(fm_k2)$cond[[1]]
sds_k2  <- attr(vc_k2, "stddev")
cors_k2 <- attr(vc_k2, "correlation")

test_that("indisting k=1: log-likelihood matches homcs", {
  expect_equal(as.numeric(logLik(fm_k1)), as.numeric(logLik(fm_k1_homcs)),
               tolerance = 1e-4)
})

test_that("indisting k=2 (4x4): symmetry constraints satisfied", {
  expect_equal(unname(sds_k2[1]),    unname(sds_k2[2]),    tolerance = 1e-4)
  expect_equal(unname(sds_k2[3]),    unname(sds_k2[4]),    tolerance = 1e-4)
  expect_equal(unname(cors_k2[1,3]), unname(cors_k2[2,4]), tolerance = 1e-4)
  expect_equal(unname(cors_k2[1,4]), unname(cors_k2[2,3]), tolerance = 1e-4)
})

test_that("indisting k=2: theta count equals k*(k+1) + 2 homcs = 8", {
  expect_equal(length(glmmTMB::getME(fm_k2, "theta")), 8L)
})

test_that("indisting: blockCode is 16", {
  vc <- getFromNamespace(".valid_covstruct", "glmmTMB")
  expect_equal(unname(vc["indisting"]), 16L)
})

test_that("indisting: non-factor member variable gives informative error", {
  sleepstudy2$member_num <- as.integer(sleepstudy2$member)
  expect_error(
    glmmTMB(Reaction ~ Days +
              indisting(0 + member_num + member_num:Days | dyad),
            data = sleepstudy2, REML = TRUE),
    "must be a factor variable"
  )
})

test_that("indisting: missing intercept suppression gives informative error", {
  expect_error(
    glmmTMB(Reaction ~ Days +
              indisting(member + member:Days | dyad),
            data = sleepstudy2, REML = TRUE),
    "must not include an intercept"
  )
})
