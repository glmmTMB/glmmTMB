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
