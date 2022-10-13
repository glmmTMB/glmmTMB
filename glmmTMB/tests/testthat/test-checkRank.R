stopifnot(require("testthat"),
          require("glmmTMB"))

set.seed(22380)
dat <- data.frame(
	y = rnorm(100),
	x1 = rnorm(100),
	x2 = rnorm(100)
)
# create two linearly dependent predictor variables
dat$x3 <- with(dat, x1 + 3*x2)
dat$x4 <- with(dat, -1*x1 + x2 + 5*x3)

test_that("error messages for non-identifiable fixed effects", {
    expect_error(
        glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='stop')),
        "fixed effects in conditional model are rank deficient"
    )
    expect_error(
        glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='stop')),
        "fixed effects in zero-inflation model are rank deficient"
    )
    expect_error(
        glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='stop')),
        "fixed effects in dispersion model are rank deficient"
    )
})

test_that("warning messages for non-identifiable fixed effects", {
    expect_warning(
        glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warning')),
        "fixed effects in conditional model are rank deficient"
    )
    expect_warning(
        glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warning')),
        "fixed effects in zero-inflation model are rank deficient"
    )
    expect_warning(
        glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warning')),
        "fixed effects in dispersion model are rank deficient"
    )
})

test_that("messages messages for non-identifiable fixed effects", {
    expect_message(
        m1 <- glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')),
        "dropping columns.*conditional")
    expect_equal(length(fixef(m1)$cond), 3L)
    expect_message(
        m1 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')),
        "dropping columns.*zero-inflation")
    expect_equal(length(fixef(m1)$zi), 3L)
    expect_message(
        m1 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')),
        "dropping columns.*dispersion")
    expect_equal(length(fixef(m1)$disp), 3L)
})

# FIXME: what tests can be performed for rank_check="skip"? maybe just check that 'skip' and 'warn' give equivalent results?

