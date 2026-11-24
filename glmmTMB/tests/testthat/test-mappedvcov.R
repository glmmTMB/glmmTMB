## downstream methods (emmeans, car::Anova) for models with
## fixed-effect coefficients fixed via 'map'

stopifnot(require("testthat"),
          require("glmmTMB"))

data(Salamanders, package = "glmmTMB")

## fix the intercept to its (rounded) unmapped estimate; spp/mined free
fit_map <- glmmTMB(count ~ mined + spp, family = poisson,
                   data = Salamanders,
                   start = list(beta = c(-1, rep(0, 7))),
                   map = list(beta = factor(c(NA, 1:7))))
fit_free <- glmmTMB(count ~ mined + spp, family = poisson,
                    data = Salamanders)

test_that("pad_mapped_vcov aligns vcov with fixef for mapped models", {
    Vred <- vcov(fit_map, include_nonest = FALSE)$cond
    Vpad <- glmmTMB:::pad_mapped_vcov(fit_map, Vred, "cond")
    expect_identical(dim(Vpad), rep(length(fixef(fit_map)$cond), 2L))
    expect_equal(unname(Vpad["(Intercept)", ]),
                 rep(0, ncol(Vpad)))
    ## full-size input: NA rows for the mapped coefficient become zeros
    Vfull <- vcov(fit_map, include_nonest = TRUE)$cond
    Vpad2 <- glmmTMB:::pad_mapped_vcov(fit_map, Vfull, "cond")
    expect_false(anyNA(Vpad2))
    expect_equal(Vpad, Vpad2)
    ## no-op for unmapped models
    Vfree <- vcov(fit_free)$cond
    expect_identical(glmmTMB:::pad_mapped_vcov(fit_free, Vfree, "cond"),
                     Vfree)
})

test_that("pad_mapped_vcov lifts the reduced vcov for coefficients tied via 'map' (GH #1340)", {
    ## two spp coefficients tied: 8 nominal coefficients, 7 estimated
    fit_tie <- glmmTMB(count ~ mined + spp, family = poisson,
                       data = Salamanders,
                       map = list(beta = factor(c(1, 2, 3, 3, 4, 5, 6, 7))))
    Vred <- vcov(fit_tie, include_nonest = FALSE)$cond
    expect_identical(dim(Vred), c(7L, 7L))
    Vpad <- glmmTMB:::pad_mapped_vcov(fit_tie, Vred, "cond")
    expect_identical(dim(Vpad), c(8L, 8L))
    expect_identical(rownames(Vpad), names(fixef(fit_tie)$cond))
    ## the tied pair shares one parameter: identical rows, correlation 1 ...
    expect_equal(Vpad[3, ], Vpad[4, ], ignore_attr = TRUE)
    ## ... and the result is what vcov(include_nonest = TRUE) reports
    expect_equal(Vpad, as.matrix(vcov(fit_tie)$cond), ignore_attr = TRUE)
    skip_if_not_installed("emmeans")
    ## used to fail in emmeans:::.qf.non0 ("logical subscript too long")
    em <- summary(emmeans::emmeans(fit_tie, ~ spp))
    expect_false(anyNA(em$SE))
    ## the tied pair (sppPR, sppDM) gives identical standard errors
    expect_equal(em$SE[em$spp == "PR"], em$SE[em$spp == "DM"])
})

test_that("emmeans works with a mapped fixed-effect coefficient", {
    skip_if_not_installed("emmeans")
    em <- summary(emmeans::emmeans(fit_map, ~ mined))
    expect_false(anyNA(em$SE))
    ## contrasts only involve estimated coefficients, so they should be
    ## close to the freely-estimated model's contrasts
    ec <- summary(emmeans::contrast(emmeans::emmeans(fit_map, ~ mined),
                                    "pairwise"))
    expect_equal(ec$estimate,
                 -unname(fixef(fit_map)$cond["minedno"]),
                 tolerance = 1e-6)
})

test_that("car::Anova works with a mapped fixed-effect coefficient", {
    skip_if_not_installed("car")
    ## type II: no NA/singular failures for free terms
    a2 <- car::Anova(fit_map)
    expect_false(anyNA(a2[["Chisq"]]))
    ## 1-df term consistent with squared z statistic
    z <- summary(fit_map)$coefficients$cond["minedno", "z value"]
    expect_equal(a2["mined", "Chisq"], z^2, tolerance = 1e-6)
    ## type III: mapped-intercept-only hypotheses yield NA, others finite
    a3 <- car::Anova(fit_map, type = 3)
    expect_true(is.na(a3["(Intercept)", "Chisq"]))
    expect_false(anyNA(a3[c("mined", "spp"), "Chisq"]))
    ## a user-supplied vcov. must be used verbatim (not padded/zeroed):
    ## with the reduced matrix supplied explicitly, Anova sees a
    ## dimension mismatch rather than silently rewriting it
    Vred <- vcov(fit_map, include_nonest = FALSE)$cond
    expect_error(suppressWarnings(car::Anova(fit_map, vcov. = Vred)))
})

test_that("car::Anova tolerates a user-supplied vcov with NA variances", {
    skip_if_not_installed("car")
    ## a full-size vcov (include_nonest = TRUE) keeps NA rows/cols for the
    ## mapped coefficient. Supplied verbatim, diag(vcov.) == 0 is NA, which
    ## previously made any(zv) return NA and errored the type II/III filters
    ## ("missing value where TRUE/FALSE needed"). It should now run and,
    ## since the matrix is used verbatim, propagate NA rather than crashing.
    Vfull <- vcov(fit_map, include_nonest = TRUE)$cond
    expect_true(anyNA(diag(Vfull)))
    a2 <- expect_no_error(car::Anova(fit_map, vcov. = Vfull))
    a3 <- expect_no_error(car::Anova(fit_map, type = 3, vcov. = Vfull))
    expect_true(all(is.na(a2[["Chisq"]])))
    expect_true(all(is.na(a3[["Chisq"]])))
})

test_that("vcov(full = TRUE) handles a mapped family parameter (psi)", {
    ## the psi block was edited from the wrong component list entry, so a
    ## mapped psi left the name list unchanged and the full vcov either
    ## silently dropped the psi row or failed to expand
    set.seed(101)
    dd <- data.frame(x = rnorm(200))
    dd$y <- 1 + 0.5 * dd$x + rt(200, df = 5)
    fit_t <- glmmTMB(y ~ x, data = dd, family = t_family(),
                     map = list(psi = factor(NA)), start = list(psi = log(5)))
    V <- vcov(fit_t, full = TRUE)
    expect_identical(unname(rownames(V)),
                     c("(Intercept)", "x", "disp~(Intercept)",
                       "Student-t df"))
    expect_true(all(is.na(V["Student-t df", ])))
    expect_false(anyNA(V[1:3, 1:3]))
    Vest <- vcov(fit_t, full = TRUE, include_nonest = FALSE)
    expect_identical(unname(rownames(Vest)), unname(rownames(V))[1:3])
    expect_equal(unname(Vest), unname(V[1:3, 1:3]))
})

test_that("vcov(full = TRUE) handles a mapped dispersion random effect", {
    ## the thetadisp map used to hit an NA index in par_components and
    ## make vcov() and summary() fail; the cond and disp theta rows share a name, so
    ## the full vcov must be filled by position
    skip_if_not_installed("lme4")
    data("sleepstudy", package = "lme4")
    fit_d <- glmmTMB(Reaction ~ Days + (1 | Subject),
                     dispformula = ~ 1 + (1 | Subject),
                     data = sleepstudy,
                     map = list(thetadisp = factor(NA)),
                     start = list(thetadisp = 0))
    V <- as.matrix(vcov(fit_d, full = TRUE))
    expect_identical(unname(rownames(V)),
                     c("(Intercept)", "Days", "disp~(Intercept)",
                       "theta_1|Subject.1", "theta_1|Subject.1"))
    expect_false(anyNA(V[1:4, 1:4]))
    expect_true(all(is.na(V[5, ])))
})

test_that("mapping a non-intercept coefficient also works", {
    skip_if_not_installed("car")
    fit_map2 <- glmmTMB(count ~ mined + spp, family = poisson,
                        data = Salamanders,
                        start = list(beta = c(0, 1, rep(0, 6))),
                        map = list(beta = factor(c(1, NA, 2:7))))
    ## the term containing only the mapped coefficient is untestable -> NA
    a2 <- car::Anova(fit_map2)
    expect_true(is.na(a2["mined", "Chisq"]))
    expect_false(is.na(a2["spp", "Chisq"]))
})

test_that("Wald confint() refuses a family parameter by name with random effects", {
    ## the Wald rows put the family parameter before the random-effect
    ## SDs, so a name resolved to its vcov(full = TRUE) position would
    ## select a random-effect row
    set.seed(5)
    d <- data.frame(g = factor(rep(1:10, each = 20)), x = rnorm(200))
    d$y <- 1 + d$x + rep(rnorm(10), each = 20) + rt(200, df = 5)
    ft <- glmmTMB(y ~ x + (1 | g), data = d, family = t_family(),
                  map = list(psi = factor(NA)), start = list(psi = log(5)))
    expect_error(confint(ft, parm = "Student-t df", include_nonest = TRUE),
                 "not supported for this family")
    ## other names still select their rows
    ci <- confint(ft)
    expect_equal(confint(ft, parm = "x"), ci["x", , drop = FALSE])
})

## 'map' factors whose levels are not in coefficient order. TMB orders the
## estimated parameters by factor level of the map, not by coefficient;
## vcov() used to label that order with the coefficient names, so
## summary() showed swapped standard errors. The maps are cycles rather
## than swaps, so that a map read in the wrong direction fails as well;
## each map is compared with a fit without 'map'
data(sleepstudy, package = "lme4")
ss <- sleepstudy
ss$Days2 <- ss$Days^2
ss$Days3 <- ss$Days^3
ss$Dcopy <- ss$Days
ss$g <- factor(rep(c("a", "b", "c"), length.out = nrow(ss)))
m0 <- glmmTMB(Reaction ~ Days + Days2 + (1 | Subject), data = ss, REML = TRUE)
## the same model, with the parameters relabelled
m1 <- update(m0, map = list(beta = factor(c(3, 1, 2))))
## a permutation together with a fixed coefficient: Days2 fixed at 0 is
## the model without Days2, and Days3 fixed at 0 is m0
m2 <- update(m0, map = list(beta = factor(c(2, 1, NA))),
             start = list(beta = c(0, 0, 0)))
m3 <- glmmTMB(Reaction ~ Days + (1 | Subject), data = ss, REML = TRUE)
m4 <- glmmTMB(Reaction ~ Days + Days2 + Days3 + (1 | Subject), data = ss,
              REML = TRUE, map = list(beta = factor(c(3, 1, 2, NA))),
              start = list(beta = c(0, 0, 0, 0)))
## the same in a rank-deficient model (Dcopy is dropped, and the map
## refers to the remaining columns)
r1 <- suppressMessages(glmmTMB(Reaction ~ Days + Dcopy + Days2 + (1 | Subject),
                               data = ss, REML = TRUE,
                               map = list(beta = factor(c(3, 1, 2)))))
r2 <- suppressMessages(glmmTMB(Reaction ~ Days + Dcopy + Days2 + (1 | Subject),
                               data = ss, REML = TRUE,
                               map = list(beta = factor(c(2, 1, NA))),
                               start = list(beta = c(0, 0, 0))))

test_that("vcov() follows the coefficients when 'map' reorders the estimated parameters", {
    expect_equal(c(logLik(m1)), c(logLik(m0)), tolerance = 1e-8)
    expect_equal(vcov(m1)$cond, vcov(m0)$cond, tolerance = 1e-6)
    expect_equal(coef(summary(m1))$cond[, "Std. Error"],
                 coef(summary(m0))$cond[, "Std. Error"], tolerance = 1e-6)
    ## the reduced matrix stays in map-level order; pad_mapped_vcov()
    ## lifts it to the coefficients
    Vred <- vcov(m1, include_nonest = FALSE)$cond
    expect_equal(glmmTMB:::pad_mapped_vcov(m1, Vred, "cond"),
                 vcov(m0)$cond, tolerance = 1e-6, ignore_attr = TRUE)
    V2 <- vcov(m2)$cond
    expect_true(all(is.na(V2["Days2", ])) && all(is.na(V2[, "Days2"])))
    expect_equal(V2[1:2, 1:2], vcov(m3)$cond, tolerance = 1e-6)
    V4 <- vcov(m4)$cond
    expect_true(all(is.na(V4["Days3", ])) && all(is.na(V4[, "Days3"])))
    expect_equal(V4[1:3, 1:3], vcov(m0)$cond, tolerance = 1e-6)
})

test_that("vcov() follows the coefficients in a rank-deficient model, without a warning", {
    expect_no_warning(vcov(r1))
    Vr1 <- suppressWarnings(vcov(r1)$cond)
    expect_true(all(is.na(Vr1["Dcopy", ])) && all(is.na(Vr1[, "Dcopy"])))
    expect_equal(Vr1[-3, -3], vcov(m0)$cond, tolerance = 1e-6)
    expect_no_warning(vcov(r2))
    Vr2 <- suppressWarnings(vcov(r2)$cond)
    expect_true(all(is.na(Vr2[c("Dcopy", "Days2"), ])))
    expect_equal(Vr2[1:2, 1:2], vcov(m3)$cond, tolerance = 1e-6)
    ## ... nor for a rank-deficient dispersion model without any map
    rd <- suppressMessages(glmmTMB(Reaction ~ Days + (1 | Subject),
                                   dispformula = ~ Days + Dcopy,
                                   data = ss, REML = TRUE))
    expect_no_warning(vcov(rd))
})

test_that("vcov() follows the coefficients in the dispersion and zero-inflation models", {
    d0 <- glmmTMB(Reaction ~ Days + (1 | Subject), dispformula = ~ g,
                  data = ss, REML = TRUE)
    d1 <- update(d0, map = list(betadisp = factor(c(3, 1, 2))))
    expect_equal(vcov(d1)$disp, vcov(d0)$disp, tolerance = 1e-6)
    set.seed(3)
    ss$cnt <- rpois(nrow(ss), 2)
    ss$cnt[sample(nrow(ss), 40)] <- 0
    z0 <- glmmTMB(cnt ~ Days + (1 | Subject), ziformula = ~ g,
                  family = poisson, data = ss)
    z1 <- update(z0, map = list(betazi = factor(c(3, 1, 2))))
    expect_equal(vcov(z1)$zi, vcov(z0)$zi, tolerance = 1e-6)
})

test_that("bread() and estfun() name the estimated parameters after the coefficients", {
    ## the estimated parameters in the order TMB uses (here the
    ## coefficients 2, 3, 1)
    b0 <- update(m0, REML = FALSE)
    b1 <- update(b0, map = list(beta = factor(c(3, 1, 2))))
    o <- order(c(3, 1, 2))
    hs <- function(m) bread(m) %*% crossprod(estfun(m)) %*% bread(m)
    expect_equal(unname(hs(b1)), unname(hs(b0)[o, o]), tolerance = 1e-6)
    expect_equal(bread(b1), bread(b0)[o, o], tolerance = 1e-6)
    expect_equal(estfun(b1), estfun(b0)[, o], tolerance = 1e-6)
    ## a tie out of order: Days and Days2 share the coefficient of DD
    bt <- update(b0, map = list(beta = factor(c(2, 1, 1))))
    ss$DD <- ss$Days + ss$Days2
    btw <- glmmTMB(Reaction ~ DD + (1 | Subject), data = ss, REML = FALSE)
    expect_identical(rownames(bread(bt)), c("Days", "(Intercept)"))
    expect_equal(unname(bread(bt)[c("(Intercept)", "Days"), c("(Intercept)", "Days")]),
                 unname(bread(btw)), tolerance = 1e-6)
    expect_identical(colnames(estfun(bt)), rownames(bread(bt)))
    expect_equal(unname(hs(bt)[c("(Intercept)", "Days"), c("(Intercept)", "Days")]),
                 unname(hs(btw)), tolerance = 1e-6)
    ## a fixed coefficient: Days2 fixed at 0 is the model without Days2
    bp <- update(b0, map = list(beta = factor(c(2, 1, NA))),
                 start = list(beta = c(0, 0, 0)))
    bw <- update(b0, . ~ . - Days2)
    expect_identical(colnames(estfun(bp)), rownames(bread(bp)))
    expect_equal(estfun(bp)[, c("(Intercept)", "Days")], estfun(bw), tolerance = 1e-6)
    ## no rows for coefficients dropped for rank deficiency (a column in
    ## the middle) ...
    bn <- suppressMessages(glmmTMB(Reaction ~ Days + Dcopy + Days2 + (1 | Subject),
                                   data = ss, REML = FALSE))
    expect_equal(bread(bn), bread(b0), tolerance = 1e-6)
    expect_equal(estfun(bn), estfun(b0), tolerance = 1e-6)
    ## ... or fixed via 'map'
    ba <- update(b0, map = list(beta = factor(c(NA, NA, NA))),
                 start = list(beta = unname(fixef(b0)$cond)))
    expect_identical(dim(bread(ba)), c(0L, 0L))
    expect_identical(dim(bread(ba, rawnames = TRUE)), c(0L, 0L))
    expect_identical(ncol(estfun(ba)), 0L)
    ## ordinal: the intercept is fixed via an internal 'map'
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    bo <- glmmTMB(rating ~ temp + contact + (1 | judge), data = wine,
                  family = ordinal())
    expect_identical(rownames(bread(bo)), c("tempwarm", "contactyes"))
    expect_identical(colnames(estfun(bo)), rownames(bread(bo)))
})

test_that("sandwich(), vcovHC() and vcov(sandwich = TRUE) follow the coefficients with a 'map'", {
    ## bread(full = TRUE) and estfun(full = TRUE) work in the space of the
    ## estimated parameters, so the sandwich exists whenever estfun()
    ## does; each fit is compared with a fit without 'map' (same model)
    b0 <- update(m0, REML = FALSE)
    sw <- function(m) list(sandwich = glmmTMB::sandwich(m),
                           vcovHC = sandwich::vcovHC(m),
                           vcov = vcov(m, sandwich = TRUE)$cond)
    cmp <- function(m, twin, nm, tnm = nm) {
        s <- sw(m); st <- sw(twin)
        for (k in names(s)) {
            expect_equal(unname(s[[k]][nm, nm]), unname(st[[k]][tnm, tnm]),
                         tolerance = 1e-6, info = k)
        }
    }
    nm <- c("(Intercept)", "Days", "Days2")
    ## a permutation
    cmp(update(b0, map = list(beta = factor(c(3, 1, 2)))), b0, nm)
    ## a tie out of order: Days and Days2 share the coefficient of DD
    ss$DD <- ss$Days + ss$Days2
    bt <- update(b0, map = list(beta = factor(c(2, 1, 1))))
    btw <- glmmTMB(Reaction ~ DD + (1 | Subject), data = ss, REML = FALSE)
    cmp(bt, btw, c("(Intercept)", "Days"), c("(Intercept)", "DD"))
    ## the parameters other than beta in the full sandwich
    expect_equal(unname(diag(glmmTMB::sandwich(bt, full = TRUE))[-(1:2)]),
                 unname(diag(glmmTMB::sandwich(btw, full = TRUE))[-(1:2)]),
                 tolerance = 1e-6)
    ## a fixed coefficient: Days2 fixed at 0 is the model without Days2
    bp <- update(b0, map = list(beta = factor(c(2, 1, NA))),
                 start = list(beta = c(0, 0, 0)))
    cmp(bp, update(b0, . ~ . - Days2), c("(Intercept)", "Days"))
    ## a coefficient dropped for rank deficiency
    bn <- suppressMessages(glmmTMB(Reaction ~ Days + Dcopy + Days2 + (1 | Subject),
                                   data = ss, REML = FALSE))
    cmp(bn, b0, nm)
    ## full = TRUE: the zi and disp blocks, theta and psi are named in the
    ## same way, and estfun() uses the names of bread()
    d0 <- glmmTMB(Reaction ~ Days + (1 | Subject), dispformula = ~ g,
                  data = ss, REML = FALSE)
    d1 <- update(d0, map = list(betadisp = factor(c(3, 1, 2))))
    nd <- c("disp~(Intercept)", "disp~gb", "disp~gc")
    expect_equal(glmmTMB::sandwich(d1, full = TRUE)[nd, nd],
                 glmmTMB::sandwich(d0, full = TRUE)[nd, nd], tolerance = 1e-6)
    expect_identical(colnames(estfun(d1, full = TRUE)), rownames(bread(d1, full = TRUE)))
    expect_identical(colnames(estfun(bt, full = TRUE)), rownames(bread(bt, full = TRUE)))
    ## a tie of standard deviations out of order: both slopes share the
    ## parameter of map level 1
    f3 <- Reaction ~ Days + (1 | Subject) + (0 + Days | Subject) + (0 + Days2 | Subject)
    tt <- glmmTMB(f3, data = ss, REML = FALSE,
                  map = list(theta = factor(c(2, 1, 1))))
    ttw <- update(tt, map = list(theta = factor(c(1, 2, 2))))
    nt <- c("theta_1|Subject.1", "theta_0+Days|Subject.1")
    expect_equal(glmmTMB::sandwich(tt, full = TRUE)[nt, nt],
                 glmmTMB::sandwich(ttw, full = TRUE)[nt, nt], tolerance = 1e-6)
    expect_identical(colnames(estfun(tt, full = TRUE)), rownames(bread(tt, full = TRUE)))
    ## without 'map' the names are those of getParnames(), except that a
    ## dispformula without fixed effects has no "disp~" row
    expect_identical(rownames(bread(d0, full = TRUE)),
                     unlist(glmmTMB:::getParnames(d0, full = TRUE), use.names = FALSE))
    ss$off <- log(30)
    dz <- glmmTMB(Reaction ~ Days + (1 | Subject), dispformula = ~ 0 + offset(off),
                  data = ss, REML = FALSE)
    expect_identical(rownames(bread(dz, full = TRUE)),
                     c("(Intercept)", "Days", "theta_1|Subject.1"))
})

test_that("emmeans() follows the coefficients when 'map' reorders them", {
    skip_if_not_installed("emmeans")
    at <- list(Days = c(0, 7), Days2 = c(0, 49))
    expect_equal(summary(emmeans::emmeans(m1, ~ Days + Days2, at = at))$SE,
                 summary(emmeans::emmeans(m0, ~ Days + Days2, at = at))$SE,
                 tolerance = 1e-6)
})

test_that("car::Anova() follows the coefficients when 'map' reorders them", {
    skip_if_not_installed("car")
    expect_equal(car::Anova(m1)$Chisq, car::Anova(m0)$Chisq, tolerance = 1e-6)
    expect_equal(car::Anova(m2)["Days", "Chisq"],
                 car::Anova(m3)["Days", "Chisq"], tolerance = 1e-6)
    expect_equal(car::Anova(r1)[c("Days", "Days2"), "Chisq"],
                 car::Anova(m0)$Chisq, tolerance = 1e-6)
    expect_equal(car::Anova(r2)["Days", "Chisq"],
                 car::Anova(m3)["Days", "Chisq"], tolerance = 1e-6)
})
