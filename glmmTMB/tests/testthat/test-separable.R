stopifnot(require("testthat"),
          require("glmmTMB"))

make_sep_dat <- function(n_member = 2, n_time = 3, n_group = 2, reps = FALSE) {
    args <- list(member = factor(paste0("m", seq_len(n_member))),
                 time = factor(seq_len(n_time)),
                 group = factor(seq_len(n_group)))
    if (reps) args$rep <- 1:2
    dd <- do.call(expand.grid, args)
    dd$y <- seq_len(nrow(dd)) / nrow(dd)
    dd
}

ar1_to_theta <- function(phi) phi / sqrt(1 - phi^2)

homcs_to_theta <- function(rho, n_member) {
    lower <- -1 / (n_member - 1)
    qlogis((rho - lower) / (1 - lower))
}

joint_nll_at <- function(form, dd, theta, b, sigma = 2) {
    fit <- glmmTMB(form, data = dd,
                   start = list(beta = 0, betadisp = log(sigma),
                                theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))
    p <- fit$obj$env$last.par.best
    p[names(p) == "beta"] <- 0
    p[names(p) == "betadisp"] <- log(sigma)
    p[names(p) == "b"] <- b
    p[names(p) == "theta"] <- theta
    fit$obj$env$f(p)
}

fit_fixed_theta <- function(form, dd, theta) {
    glmmTMB(form, data = dd,
            start = list(theta = theta),
            map = list(theta = factor(rep(NA, length(theta)))))
}

make_sep_case <- function(struc = c("homcs", "us"), reversed = FALSE,
                          n_member = 2, n_time = 3) {
    struc <- match.arg(struc)
    rho <- if (struc == "homcs") 0.3 else -0.25
    phi <- if (struc == "homcs") 0.4 else 0.5
    sd <- if (struc == "homcs") 2 else seq(0.8, 1.2, length.out = n_member)

    R_member <- if (struc == "homcs") {
        M <- matrix(rho, n_member, n_member)
        diag(M) <- 1
        M
    } else {
        outer(seq_len(n_member), seq_len(n_member),
              function(i, j) ifelse(i == j, 1, 0.35^abs(i - j)))
    }
    R_time <- outer(seq_len(n_time), seq_len(n_time),
                    function(i, j) phi^abs(i - j))
    member_theta <- switch(struc,
        homcs = c(log(sd), homcs_to_theta(rho, n_member)),
        us = c(log(sd), put_cor(R_member))
    )

    if (reversed) {
        form <- switch(struc,
            homcs = y ~ 1 + separable(ar1(0 + time) %x% homcs(0 + member) | group),
            us = y ~ 1 + separable(ar1(0 + time) %x% us(0 + member) | group)
        )
        dense_form <- y ~ 1 + us(sepgrid(time, member) + 0 | group)
        theta <- c(ar1_to_theta(phi), member_theta)
        R_full <- kronecker(R_member, R_time)
        sd_full <- if (length(sd) == 1) rep(sd, n_member * n_time) else rep(sd, each = n_time)
        codes <- unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[[struc]]))
        kinds <- c(2L, 1L)
        scale_spec <- 1L
    } else {
        form <- switch(struc,
            homcs = y ~ 1 + separable(homcs(0 + member) %x% ar1(0 + time) | group),
            us = y ~ 1 + separable(us(0 + member) %x% ar1(0 + time) | group)
        )
        dense_form <- y ~ 1 + us(sepgrid(member, time) + 0 | group)
        theta <- c(member_theta, ar1_to_theta(phi))
        R_full <- kronecker(R_time, R_member)
        sd_full <- if (length(sd) == 1) rep(sd, n_member * n_time) else rep(sd, n_time)
        codes <- unname(c(.valid_covstruct[[struc]], .valid_covstruct[["ar1"]]))
        kinds <- c(1L, 2L)
        scale_spec <- 0L
    }

    list(form = form, dense_form = dense_form, theta = theta,
         theta_dense = c(log(sd_full), put_cor(R_full)),
         R_full = R_full, sd_full = sd_full,
         codes = codes, kinds = kinds, scale_spec = scale_spec,
         n_member = n_member, n_time = n_time)
}

expect_separable_vc <- function(case) {
    dd <- make_sep_dat(n_member = case$n_member, n_time = case$n_time, reps = TRUE)
    fit <- fit_fixed_theta(case$form, dd, case$theta)
    restruc <- fit$modelInfo$reStruc$condReStruc[[1]]
    vc <- VarCorr(fit)$cond[[1]]

    expect_equal(restruc$sepCodes, case$codes)
    expect_equal(restruc$sepDensityKinds, case$kinds)
    expect_equal(restruc$sepDispatch, 1L)
    expect_equal(restruc$sepScaleMode, 1L)
    expect_equal(restruc$sepScaleSpec, case$scale_spec)
    expect_equal(unname(attr(vc, "stddev")), case$sd_full, tolerance = 1e-6)
    expect_equal(unname(attr(vc, "correlation")), case$R_full, tolerance = 1e-6)
}

expect_separable_dense_nll <- function(case) {
    dd <- make_sep_dat(n_member = case$n_member, n_time = case$n_time, n_group = 1)
    dd$y <- 0
    b <- seq(-0.4, 0.5, length.out = length(case$sd_full))

    sep_nll <- joint_nll_at(case$form, dd, case$theta, b)
    dense_nll <- joint_nll_at(case$dense_form, dd, case$theta_dense, b)

    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)
}

test_that("sepgrid builds complete two-dimensional levels", {
    member <- factor(c("A", "B"), levels = c("A", "B"))
    time <- factor(c(1, 3), levels = 1:3)
    grid <- sepgrid(member, time)

    expect_equal(nlevels(grid), 6)
    expect_equal(unname(parseNumLevels(levels(grid))),
                 unname(as.matrix(expand.grid(1:2, 1:3))),
                 check.attributes = FALSE)
})

test_that("glmmTMB preserves unused sepgrid levels by default", {
    dd <- expand.grid(member = factor(c("A", "B")),
                      time = factor(c(1, 3), levels = 1:3),
                      group = factor(1:2),
                      rep = 1:2)
    dd$fixed_factor <- factor(rep(c("a", "b"), length.out = nrow(dd)),
                              levels = c("a", "b", "unused"))
    dd$y <- 0

    fit <- glmmTMB(y ~ fixed_factor +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepDims, c(2L, 3L))
    expect_equal(levels(fit$fr[["sepgrid(member, time)"]]),
                 c("(1,1)", "(2,1)", "(1,2)", "(2,2)", "(1,3)", "(2,3)"))
    expect_equal(levels(fit$fr$fixed_factor), c("a", "b"))
})

test_that("separable metadata handles product order and un alias", {
    dd <- make_sep_dat()

    h <- glmmTMB(y ~ 1 +
                     separable(ar1(0 + time) %x% homcs(0 + member) | group),
                 data = dd, doFit = FALSE)
    u <- glmmTMB(y ~ 1 +
                     separable(un(0 + member) %x% ar1(0 + time) | group),
                 data = dd, doFit = FALSE)

    expect_equal(unname(h$condReStruc[[1]]$blockCode),
                 unname(.valid_covstruct[["separable"]]))
    expect_equal(h$condReStruc[[1]]$blockNumTheta, 3)
    expect_equal(h$condReStruc[[1]]$sepDims, c(3L, 2L))
    expect_equal(h$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[["homcs"]])))
    expect_equal(h$condReStruc[[1]]$sepDensityKinds, c(2L, 1L))
    expect_equal(h$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleSpec, 1L)

    expect_equal(u$condReStruc[[1]]$blockNumTheta, 4)
    expect_equal(u$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["us"]], .valid_covstruct[["ar1"]])))
    expect_equal(u$condReStruc[[1]]$sepDensityKinds, c(1L, 2L))
    expect_equal(u$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(u$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(u$condReStruc[[1]]$sepScaleSpec, 0L)
})

test_that("separable parser keeps structured metadata in splitForm payload", {
    f <- y ~ 1 +
        separable(homcs(0 + member) %x% ar1(0 + time) | group,
                  scale = homcs(0 + member))
    g <- glmmTMB:::rewrite_separable_formula(f)
    ss <- reformulas::splitForm(g, specials = c(names(.valid_covstruct), "s"))
    spec <- eval(ss$reTrmAddArgs[[1]][[2]])

    expect_null(attr(g, "separable_specs"))
    expect_equal(spec$grid, c("member", "time"))
    expect_equal(unname(spec$margins[, "struc"]), c("homcs", "ar1"))
    expect_equal(unname(spec$scale[, "struc"]), "homcs")
    expect_equal(unname(spec$scale[, "var"]), "member")
    expect_equal(ss$reTrmClasses, "separable")
    expect_equal(deparse(ss$reTrmFormulas[[1]]),
                 "sepgrid(member, time) + 0 | group")
    expect_equal(length(ss$reTrmAddArgs[[1]]), 2)
})

test_that("separable product syntax rejects slope-like margins for now", {
    dd <- make_sep_dat()
    dd$x <- seq_len(nrow(dd))

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member + x) %x% ar1(0 + time) | group,
                              scale = us(0 + member + x)),
                data = dd, doFit = FALSE),
        "exactly one no-intercept variable"
    )
})

test_that("separable supports explicit scale margin selection", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1 +
                       separable(ar1(0 + time) %x% un(0 + member) | group,
                                 scale = un(0 + member)),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[["us"]])))
    expect_equal(fit$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(fit$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(fit$condReStruc[[1]]$sepScaleSpec, 1L)
    expect_equal(fit$condReStruc[[1]]$blockNumTheta, 4)
})

test_that("multiple separable terms keep their metadata order", {
    dd <- expand.grid(member = factor(c("A", "B")),
                      time = factor(1:3),
                      group1 = factor(1:2),
                      group2 = factor(1:2))
    dd$y <- 0

    fit <- glmmTMB(y ~ 1 +
                       (1 | group1) +
                       separable(un(0 + member) %x% ar1(0 + time) | group1) +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group2),
                   data = dd, doFit = FALSE)

    sep_terms <- vapply(fit$condReStruc, function(x) {
        identical(unname(x$blockCode), unname(.valid_covstruct[["separable"]]))
    }, logical(1))

    expect_equal(unname(which(sep_terms)), c(2L, 3L))
    expect_equal(fit$condReStruc[[2]]$blockNumTheta, 4)
    expect_equal(fit$condReStruc[[2]]$sepCodes,
                 unname(c(.valid_covstruct[["us"]], .valid_covstruct[["ar1"]])))
    expect_equal(fit$condReStruc[[3]]$blockNumTheta, 3)
    expect_equal(fit$condReStruc[[3]]$sepCodes,
                 unname(c(.valid_covstruct[["homcs"]], .valid_covstruct[["ar1"]])))
})

test_that("separable preserves the user-facing formula in the stored call", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    ftxt <- paste(deparse(fit$call$formula), collapse = " ")
    expect_match(ftxt, "homcs\\(0 \\+ member\\)")
    expect_false(grepl("data.frame", ftxt, fixed = TRUE))
})

test_that("separable preserves user-facing zi and dispersion formulas", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1,
                   ziformula = ~ separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   dispformula = ~ separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd, doFit = FALSE)

    ztxt <- paste(deparse(fit$call$ziformula), collapse = " ")
    dtxt <- paste(deparse(fit$call$dispformula), collapse = " ")
    expect_match(ztxt, "homcs\\(0 \\+ member\\)")
    expect_match(dtxt, "homcs\\(0 \\+ member\\)")
    expect_false(grepl("data.frame", ztxt, fixed = TRUE))
    expect_false(grepl("data.frame", dtxt, fixed = TRUE))
})

test_that("separable rejects unsupported margins", {
    dd <- make_sep_dat()

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% homcs(0 + time) | group),
                data = dd, doFit = FALSE),
        "specify the scale margin"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% ar1(0 + time) | group,
                              scale = ar1(0 + time)),
                data = dd, doFit = FALSE),
        "correlation-only"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% ar1(0 + time) | group,
                              scale = homcs(0 + member)),
                data = dd, doFit = FALSE),
        "scale must match"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(ar1(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "global scale"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(us(0 + member) %x% homcs(0 + time) | group,
                              scale = us(0 + member)),
                data = dd, doFit = FALSE),
        "does not yet implement"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(foo(0 + member) %x% ar1(0 + time) | group),
                data = dd, doFit = FALSE),
        "Unsupported separable\\(\\) margin: foo"
    )
})

test_that("separable reports kronecker covariance for supported dense x ar1 pairs", {
    cases <- list(
        make_sep_case("homcs", n_time = 4),
        make_sep_case("homcs", reversed = TRUE, n_time = 4),
        make_sep_case("homcs", n_member = 5, n_time = 4),
        make_sep_case("homcs", reversed = TRUE, n_member = 5, n_time = 4),
        make_sep_case("us"),
        make_sep_case("us", reversed = TRUE)
    )
    invisible(lapply(cases, expect_separable_vc))
})

test_that("separable likelihood matches dense MVN for supported dense x ar1 pairs", {
    cases <- list(
        make_sep_case("homcs"),
        make_sep_case("homcs", reversed = TRUE),
        make_sep_case("homcs", n_member = 4),
        make_sep_case("us", n_member = 3),
        make_sep_case("us"),
        make_sep_case("us", reversed = TRUE)
    )
    invisible(lapply(cases, expect_separable_dense_nll))
})

test_that("separable prediction with newdata reports current limitation", {
    dd <- make_sep_dat()
    theta <- c(log(1), qlogis((0.2 + 1) / 2), ar1_to_theta(0.3))
    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    expect_error(predict(fit, newdata = dd[1, ]),
                 "newdata is not yet implemented")
})

test_that("separable simulation reports current limitation", {
    dd <- make_sep_dat(n_time = 2, reps = TRUE)

    theta <- c(log(1), qlogis((0.2 + 1) / 2), ar1_to_theta(0.3))
    fit <- glmmTMB(y ~ 1 +
                       separable(homcs(0 + member) %x% ar1(0 + time) | group),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    expect_error(simulate(fit, nsim = 1),
                 "simulation is not yet implemented for separable covariance structures")
})
