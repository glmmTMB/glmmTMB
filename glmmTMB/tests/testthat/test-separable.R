stopifnot(require("testthat"),
          require("glmmTMB"))

make_sep_dat <- function(n_time = 3, n_group = 2, reps = FALSE) {
    args <- list(member = factor(c("A", "B")),
                 time = factor(seq_len(n_time)),
                 group = factor(seq_len(n_group)))
    if (reps) args$rep <- 1:2
    dd <- do.call(expand.grid, args)
    dd$y <- seq_len(nrow(dd)) / nrow(dd)
    dd
}

ar1_to_theta <- function(phi) phi / sqrt(1 - phi^2)

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
                       separable(sepgrid(member, time) + 0 | group,
                                 homcs(member), ar1(time)),
                   data = dd, doFit = FALSE)

    expect_equal(fit$condReStruc[[1]]$sepDims, c(2L, 3L))
    expect_equal(levels(fit$fr[["sepgrid(member, time)"]]),
                 c("(1,1)", "(2,1)", "(1,2)", "(2,2)", "(1,3)", "(2,3)"))
    expect_equal(levels(fit$fr$fixed_factor), c("a", "b"))
})

test_that("separable metadata handles margin order and un alias", {
    dd <- make_sep_dat()

    h <- glmmTMB(y ~ 1 +
                     separable(sepgrid(member, time) + 0 | group,
                               ar1(time), homcs(member)),
                 data = dd, doFit = FALSE)
    u <- glmmTMB(y ~ 1 +
                     separable(sepgrid(member, time) + 0 | group,
                               ar1(time), un(member)),
                 data = dd, doFit = FALSE)

    expect_equal(unname(h$condReStruc[[1]]$blockCode),
                 unname(.valid_covstruct[["separable"]]))
    expect_equal(h$condReStruc[[1]]$blockNumTheta, 3)
    expect_equal(h$condReStruc[[1]]$sepDims, c(2L, 3L))
    expect_equal(h$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["homcs"]], .valid_covstruct[["ar1"]])))
    expect_equal(h$condReStruc[[1]]$sepDensityKinds, c(1L, 2L))
    expect_equal(h$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(h$condReStruc[[1]]$sepScaleSpec, 0L)

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
        separable(sepgrid(member, time) + 0 | group,
                  homcs(member), ar1(time),
                  scale = homcs(member))
    g <- glmmTMB:::rewrite_separable_formula(f)
    ss <- reformulas::splitForm(g, specials = c(names(.valid_covstruct), "s"))
    spec <- eval(ss$reTrmAddArgs[[1]][[2]])

    expect_null(attr(g, "separable_specs"))
    expect_equal(spec$grid, c("member", "time"))
    expect_equal(unname(spec$margins[, "struc"]), c("homcs", "ar1"))
    expect_equal(unname(spec$scale[, "struc"]), "homcs")
    expect_equal(unname(spec$scale[, "var"]), "member")
    expect_match(spec$user_call, "separable")
    expect_equal(ss$reTrmClasses, "separable")
    expect_equal(length(ss$reTrmAddArgs[[1]]), 2)
})

test_that("separable supports explicit scale margin selection", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1 +
                       separable(sepgrid(time, member) + 0 | group,
                                 ar1(time), un(member),
                                 scale = un(member)),
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
                       separable(sepgrid(member, time) + 0 | group1,
                                 ar1(time), un(member)) +
                       separable(sepgrid(member, time) + 0 | group2,
                                 homcs(member), ar1(time)),
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
                       separable(sepgrid(member, time) + 0 | group,
                                 homcs(member), ar1(time)),
                   data = dd, doFit = FALSE)

    ftxt <- paste(deparse(fit$call$formula), collapse = " ")
    expect_match(ftxt, "homcs\\(member\\)")
    expect_false(grepl("data.frame", ftxt, fixed = TRUE))
})

test_that("separable preserves user-facing zi and dispersion formulas", {
    dd <- make_sep_dat()

    fit <- glmmTMB(y ~ 1,
                   ziformula = ~ separable(sepgrid(member, time) + 0 | group,
                                           homcs(member), ar1(time)),
                   dispformula = ~ separable(sepgrid(member, time) + 0 | group,
                                             homcs(member), ar1(time)),
                   data = dd, doFit = FALSE)

    ztxt <- paste(deparse(fit$call$ziformula), collapse = " ")
    dtxt <- paste(deparse(fit$call$dispformula), collapse = " ")
    expect_match(ztxt, "homcs\\(member\\)")
    expect_match(dtxt, "homcs\\(member\\)")
    expect_false(grepl("data.frame", ztxt, fixed = TRUE))
    expect_false(grepl("data.frame", dtxt, fixed = TRUE))
})

test_that("separable rejects unsupported margins", {
    dd <- make_sep_dat()

    expect_error(
        glmmTMB(y ~ 1 +
                    separable(sepgrid(member, time) + 0 | group,
                              us(member), homcs(time)),
                data = dd, doFit = FALSE),
        "specify the scale margin"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(sepgrid(member, time) + 0 | group,
                              us(member), ar1(time),
                              scale = ar1(time)),
                data = dd, doFit = FALSE),
        "correlation-only"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(sepgrid(member, time) + 0 | group,
                              us(member), ar1(time),
                              scale = homcs(member)),
                data = dd, doFit = FALSE),
        "scale must match"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(sepgrid(member, time) + 0 | group,
                              ar1(member), ar1(time)),
                data = dd, doFit = FALSE),
        "global scale"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(sepgrid(member, time) + 0 | group,
                              us(member), homcs(time),
                              scale = us(member)),
                data = dd, doFit = FALSE),
        "does not yet implement"
    )
    expect_error(
        glmmTMB(y ~ 1 +
                    separable(sepgrid(member, time) + 0 | group,
                              foo(member), ar1(time)),
                data = dd, doFit = FALSE),
        "Unsupported separable\\(\\) margin: foo"
    )
})

test_that("separable homcs x ar1 reports kronecker covariance", {
    dd <- make_sep_dat(n_time = 4, reps = TRUE)

    rho <- 0.3
    phi <- 0.4
    sd <- 2
    theta <- c(log(sd), qlogis((rho + 1) / 2), ar1_to_theta(phi))

    fit <- glmmTMB(y ~ 1 +
                       separable(sepgrid(member, time) + 0 | group,
                                 homcs(member), ar1(time)),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    vc <- VarCorr(fit)$cond[[1]]
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)
    R_time <- outer(1:4, 1:4, function(i, j) phi^abs(i - j))

    expect_equal(unname(attr(vc, "stddev")), rep(sd, 8), tolerance = 1e-6)
    expect_equal(unname(attr(vc, "correlation")),
                 kronecker(R_time, R_member), tolerance = 1e-6)
})

test_that("separable homcs x ar1 likelihood matches dense MVN", {
    dd <- make_sep_dat(n_time = 3, n_group = 1)
    dd$y <- 0

    rho <- 0.3
    phi <- 0.4
    sd <- 1.7
    b <- seq(-0.4, 0.5, length.out = 6)
    theta <- c(log(sd), qlogis((rho + 1) / 2), ar1_to_theta(phi))
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)
    R_time <- outer(1:3, 1:3, function(i, j) phi^abs(i - j))
    R_full <- kronecker(R_time, R_member)
    theta_dense <- c(rep(log(sd), length(b)), put_cor(R_full))

    sep_nll <- joint_nll_at(
        y ~ 1 + separable(sepgrid(member, time) + 0 | group,
                          homcs(member), ar1(time)),
        dd, theta, b)
    dense_nll <- joint_nll_at(
        y ~ 1 + us(sepgrid(member, time) + 0 | group),
        dd, theta_dense, b)

    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)
})

test_that("separable homcs x ar1 can use reversed sepgrid order", {
    dd <- make_sep_dat(n_time = 4, reps = TRUE)

    rho <- 0.3
    phi <- 0.4
    sd <- 2
    theta <- c(ar1_to_theta(phi),
               log(sd),
               qlogis((rho + 1) / 2))

    fit <- glmmTMB(y ~ 1 +
                       separable(sepgrid(time, member) + 0 | group,
                                 ar1(time), homcs(member)),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    vc <- VarCorr(fit)$cond[[1]]
    R_time <- outer(1:4, 1:4, function(i, j) phi^abs(i - j))
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)

    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[["homcs"]])))
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepDensityKinds, c(2L, 1L))
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepScaleSpec, 1L)
    expect_equal(unname(attr(vc, "stddev")), rep(sd, 8), tolerance = 1e-6)
    expect_equal(unname(attr(vc, "correlation")),
                 kronecker(R_member, R_time), tolerance = 1e-6)
})

test_that("reversed separable homcs x ar1 likelihood matches dense MVN", {
    dd <- make_sep_dat(n_time = 3, n_group = 1)
    dd$y <- 0

    rho <- 0.3
    phi <- 0.4
    sd <- 1.7
    b <- seq(-0.4, 0.5, length.out = 6)
    theta <- c(ar1_to_theta(phi), log(sd), qlogis((rho + 1) / 2))
    R_time <- outer(1:3, 1:3, function(i, j) phi^abs(i - j))
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)
    R_full <- kronecker(R_member, R_time)
    theta_dense <- c(rep(log(sd), length(b)), put_cor(R_full))

    sep_nll <- joint_nll_at(
        y ~ 1 + separable(sepgrid(time, member) + 0 | group,
                          ar1(time), homcs(member)),
        dd, theta, b)
    dense_nll <- joint_nll_at(
        y ~ 1 + us(sepgrid(time, member) + 0 | group),
        dd, theta_dense, b)

    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)
})

test_that("separable us x ar1 reports kronecker covariance", {
    dd <- make_sep_dat(reps = TRUE)

    rho <- -0.25
    phi <- 0.5
    sd <- c(1.2, 0.8)
    theta <- c(log(sd), put_cor(rho, input_val = "vec"),
               ar1_to_theta(phi))

    fit <- glmmTMB(y ~ 1 +
                       separable(sepgrid(member, time) + 0 | group,
                                 us(member), ar1(time)),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    vc <- VarCorr(fit)$cond[[1]]
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)
    R_time <- outer(1:3, 1:3, function(i, j) phi^abs(i - j))

    expect_equal(unname(attr(vc, "stddev")), rep(sd, 3), tolerance = 1e-6)
    expect_equal(unname(attr(vc, "correlation")),
                 kronecker(R_time, R_member), tolerance = 1e-6)
})

test_that("separable us x ar1 likelihood matches dense MVN", {
    dd <- make_sep_dat(n_group = 1)
    dd$y <- 0

    rho <- -0.25
    phi <- 0.5
    sd <- c(1.2, 0.8)
    b <- seq(-0.4, 0.5, length.out = 6)
    theta <- c(log(sd), put_cor(rho, input_val = "vec"),
               ar1_to_theta(phi))
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)
    R_time <- outer(1:3, 1:3, function(i, j) phi^abs(i - j))
    R_full <- kronecker(R_time, R_member)
    sd_full <- rep(sd, 3)
    theta_dense <- c(log(sd_full), put_cor(R_full))

    sep_nll <- joint_nll_at(
        y ~ 1 + separable(sepgrid(member, time) + 0 | group,
                          us(member), ar1(time)),
        dd, theta, b)
    dense_nll <- joint_nll_at(
        y ~ 1 + us(sepgrid(member, time) + 0 | group),
        dd, theta_dense, b)

    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)
})

test_that("separable us x ar1 can use reversed sepgrid order", {
    dd <- make_sep_dat(reps = TRUE)

    rho <- -0.25
    phi <- 0.5
    sd <- c(1.2, 0.8)
    theta <- c(ar1_to_theta(phi),
               log(sd),
               put_cor(rho, input_val = "vec"))

    fit <- glmmTMB(y ~ 1 +
                       separable(sepgrid(time, member) + 0 | group,
                                 ar1(time), us(member)),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    vc <- VarCorr(fit)$cond[[1]]
    R_time <- outer(1:3, 1:3, function(i, j) phi^abs(i - j))
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)

    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepCodes,
                 unname(c(.valid_covstruct[["ar1"]], .valid_covstruct[["us"]])))
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepDensityKinds, c(2L, 1L))
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepDispatch, 1L)
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepScaleMode, 1L)
    expect_equal(fit$modelInfo$reStruc$condReStruc[[1]]$sepScaleSpec, 1L)
    expect_equal(unname(attr(vc, "stddev")), rep(sd, each = 3), tolerance = 1e-6)
    expect_equal(unname(attr(vc, "correlation")),
                 kronecker(R_member, R_time), tolerance = 1e-6)
})

test_that("reversed separable us x ar1 likelihood matches dense MVN", {
    dd <- make_sep_dat(n_group = 1)
    dd$y <- 0

    rho <- -0.25
    phi <- 0.5
    sd <- c(1.2, 0.8)
    b <- seq(-0.4, 0.5, length.out = 6)
    theta <- c(ar1_to_theta(phi),
               log(sd),
               put_cor(rho, input_val = "vec"))
    R_time <- outer(1:3, 1:3, function(i, j) phi^abs(i - j))
    R_member <- matrix(c(1, rho, rho, 1), 2, 2)
    R_full <- kronecker(R_member, R_time)
    sd_full <- rep(sd, each = 3)
    theta_dense <- c(log(sd_full), put_cor(R_full))

    sep_nll <- joint_nll_at(
        y ~ 1 + separable(sepgrid(time, member) + 0 | group,
                          ar1(time), us(member)),
        dd, theta, b)
    dense_nll <- joint_nll_at(
        y ~ 1 + us(sepgrid(time, member) + 0 | group),
        dd, theta_dense, b)

    expect_equal(unname(sep_nll), unname(dense_nll), tolerance = 1e-6)
})

test_that("separable prediction with newdata reports current limitation", {
    dd <- make_sep_dat()
    theta <- c(log(1), qlogis((0.2 + 1) / 2), ar1_to_theta(0.3))
    fit <- glmmTMB(y ~ 1 +
                       separable(sepgrid(member, time) + 0 | group,
                                 homcs(member), ar1(time)),
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
                       separable(sepgrid(member, time) + 0 | group,
                                 homcs(member), ar1(time)),
                   data = dd,
                   start = list(theta = theta),
                   map = list(theta = factor(rep(NA, length(theta)))))

    expect_error(simulate(fit, nsim = 1),
                 "simulation is not yet implemented for separable covariance structures")
})
