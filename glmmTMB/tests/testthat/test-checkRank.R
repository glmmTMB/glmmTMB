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
    # X
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
    # sparse X
    expect_error(
        glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='stop'), sparseX=c(cond=TRUE)),
        "fixed effects in conditional model are rank deficient"
    )
    expect_error(
        glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='stop'), sparseX=c(zi=TRUE)),
        "fixed effects in zero-inflation model are rank deficient"
    )
    expect_error(
        glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='stop'), sparseX=c(disp=TRUE)),
        "fixed effects in dispersion model are rank deficient"
    )
})

test_that("warning messages for non-identifiable fixed effects", {
    # X
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
    # sparse X
    expect_warning(
        glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warning'), sparseX=c(cond=TRUE)),
        "fixed effects in conditional model are rank deficient"
    )
    expect_warning(
        glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warning'), sparseX=c(zi=TRUE)),
        "fixed effects in zero-inflation model are rank deficient"
    )
    expect_warning(
        glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warning'), sparseX=c(disp=TRUE)),
        "fixed effects in dispersion model are rank deficient"
    )
})

test_that("messages for non-identifiable fixed effects", {
    # X
    expect_message(
        m1 <- glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')),
        "dropping columns.*conditional")
    expect_equal(length(fixef(m1)$cond), 5L)
    expect_equal(unname(fixef(m1)$cond[c("x3", "x4")]), rep(NA_real_, 2))

    expect_message(
        m1 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')),
        "dropping columns.*zero-inflation")
    expect_equal(length(fixef(m1)$zi), 5L)
    expect_equal(unname(fixef(m1)$zi[c("x3", "x4")]), rep(NA_real_, 2))
    expect_message(
        m1 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')),
        "dropping columns.*dispersion")
    expect_equal(length(fixef(m1)$disp), 5L)
    expect_equal(unname(fixef(m1)$disp[c("x3", "x4")]), rep(NA_real_, 2))
    ## sparse X
    ## FIXME: NA-substitution not working yet for sparse X
    ## adj.X attribute missing ... ??
    expect_message(
        m1 <- glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust'), sparseX=c(cond=TRUE)),
        "dropping columns.*conditional")
    expect_equal(length(fixef(m1)$cond), 3L)
    expect_message(
        m1 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust'), sparseX=c(zi=TRUE)),
        "dropping columns.*zero-inflation")
    expect_equal(length(fixef(m1)$zi), 3L)
    expect_message(
        m1 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust'), sparseX=c(disp=TRUE)),
        "dropping columns.*dispersion")
    expect_equal(length(fixef(m1)$disp), 3L)

})


test_that("vcov for rank-deficient models", {
    m1 <-
        suppressMessages(glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')))
    expect_equal(dim(vcov(m1)$cond), c(5, 5))
    expect_equal(dim(vcov(m1, include_nonest = FALSE)$cond), c(3, 3))
})

test_that("predict for rank-deficient models", {
    m1 <-
        suppressMessages(glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='adjust')))

    expect_equal(head(predict(m1, re.form = ~ 0), 4),
                 c(-0.00552313394691104, -0.0330178636893558, -0.214547485709826, 
                   0.127729913085515), tolerance = 1e-5)

})


test_that("equivalence between 'skip' and 'warn' when confronted with identifiable and non-identifiable fixed effects", {
    # models with no identifiability issues
    # X
    m1 <- glmmTMB(y ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='skip'))
    m2 <- glmmTMB(y ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='warn'))
    expect_equal(fixef(m1), fixef(m2))
    m1 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='skip'))
    m2 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='warn'))
    expect_equal(fixef(m1), fixef(m2))
    m1 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='skip'))
    m2 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='warn'))
    expect_equal(fixef(m1), fixef(m2))
    # sparse X
    m1 <- glmmTMB(y ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='skip'), sparseX=c(cond=TRUE))
    m2 <- glmmTMB(y ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='warn'), sparseX=c(cond=TRUE))
    expect_equal(fixef(m1), fixef(m2))
    m1 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='skip'), sparseX=c(zi=TRUE))
    m2 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='warn'), sparseX=c(zi=TRUE))
    expect_equal(fixef(m1), fixef(m2))
    m1 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='skip'), sparseX=c(disp=TRUE))
    m2 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2, data=dat, control=glmmTMBControl(rank_check='warn'), sparseX=c(disp=TRUE))
    expect_equal(fixef(m1), fixef(m2))
    # models with identifiability issues
    ## X
    cc1 <- glmmTMBControl(rank_check = 'skip', conv_check = 'skip')
    ## small-eig warning
    suppressWarnings(m1 <- glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=cc1))
    expect_warning(
        m2 <- glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warn')),
        "fixed effects in conditional model are rank deficient"
    )
    expect_equal(fixef(m1), fixef(m2))
    ## small-eig warning
    suppressWarnings(
        m1 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=cc1)
    )
    expect_warning(
        m2 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warn')),
        "fixed effects in zero-inflation model are rank deficient"
    )
    expect_equal(fixef(m1), fixef(m2))
    ## suppressWarnings() since we know this will be a non-pos-def Hessian
    m1 <- suppressWarnings(glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='skip')))
    expect_warning(
        m2 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warn')),
        "fixed effects in dispersion model are rank deficient"
    )
    expect_equal(fixef(m1), fixef(m2))
    ## sparse X
    ## small-eig warning
    suppressWarnings(m1 <- glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=cc1))
    expect_warning(
        m2 <- glmmTMB(y ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warn'), sparseX=c(cond=TRUE)),
        "fixed effects in conditional model are rank deficient"
    )
    expect_equal(fixef(m1), fixef(m2))
    ## small-eig warning
    suppressWarnings(
        m1 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control = cc1, sparseX=c(zi=TRUE))
    )
    expect_warning(
        m2 <- glmmTMB(y ~ 1, ziformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warn'), sparseX=c(zi=TRUE)),
        "fixed effects in zero-inflation model are rank deficient"
    )
    expect_equal(fixef(m1), fixef(m2))
    ## suppressWarnings() since we know this will be a non-pos-def Hessian
    m1 <- suppressWarnings(glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='skip'), sparseX=c(cond=TRUE)))
    expect_warning(
        m2 <- glmmTMB(y ~ 1, dispformula = ~ x1 + x2 + x3 + x4, data=dat, control=glmmTMBControl(rank_check='warn'), sparseX=c(cond=TRUE)),
        "fixed effects in dispersion model are rank deficient"
    )
    expect_equal(fixef(m1), fixef(m2))
})
