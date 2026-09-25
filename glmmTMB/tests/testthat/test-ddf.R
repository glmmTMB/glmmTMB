library(glmmTMB)
if (requireNamespace("pbkrtest") && requireNamespace("lme4")) {
    fm1_lmer <- lme4::lmer(formula(fm1), lme4::sleepstudy)
    fm2_lmer <- lme4::lmer(formula(fm2), lme4::sleepstudy)

    fm1 <- update(fm1, REML = TRUE)
    fm2 <- update(fm2, REML = TRUE)
    
    pbkrtest_dof <- function(m) {
        vva <- pbkrtest::vcovAdj(m)
        vv0 <- vcov(m)
        p <- length(fixef(m))
        Lmat <- diag(p)
        ddf <- apply(Lmat, 1, \(L) pbkrtest::ddf_Lb(vva, L, vv0))
        return(ddf)
    }

    test_that("approx KR df match with pbkrtest (fm1)", {
        df1 <- pbkrtest_dof(fm1_lmer)
        df2 <- dof_KR(fm1)
        expect_equal(df1, unname(c(df2)), tolerance = 1e-6)
    })

    test_that("approx KR df match with pbkrtest (fm2)", {
        df1 <- pbkrtest_dof(fm2_lmer)
        df2 <- dof_KR(fm2)
        expect_equal(df1, unname(c(df2)), tolerance = 1e-6)
    })

    test_that("KR in summary", {
        expect_identical(
            coef(summary(fm1, ddf = "kenward-roger"))$cond[,"ddf"],
            c(dof_KR(fm1))
        )
    })
}


test_that("Satt in summary", {
        expect_identical(
            unname(coef(summary(fm1, ddf = "satterthwaite"))$cond[,"ddf"]),
            c(dof_satt(fm1))
        )
    })               

test_that("Satterthwaite ddf match lmerTest (hard-coded reference values)", {
    ## reference values computed once via lmerTest 3.x (not at test time, to
    ## avoid adding a test dependency on lmerTest):
    ##   library(lmerTest)
    ##   fm1 <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = TRUE)
    ##   fm2 <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = TRUE)
    ##   coef(summary(fm1))[, "df"]; coef(summary(fm2))[, "df"]
    lmerTest_ddf_fm1 <- c(`(Intercept)` = 22.8102, Days = 161.0000)
    lmerTest_ddf_fm2 <- c(`(Intercept)` = 16.99973, Days = 16.99998)

    ## refit with REML = TRUE explicitly rather than relying on the
    ## conditional pbkrtest/lme4 block above (which only reassigns fm1/fm2 to
    ## REML fits if pbkrtest happens to be installed)
    fm1_reml <- update(fm1, REML = TRUE)
    fm2_reml <- update(fm2, REML = TRUE)

    ## tolerance of 1% comfortably covers the small, expected numerical
    ## difference between glmmTMB's finite-difference Satterthwaite
    ## implementation and lmerTest's (differences were ~0.03% in testing)
    expect_equal(unname(dof_satt(fm1_reml)), unname(lmerTest_ddf_fm1), tolerance = 1e-2)
    expect_equal(unname(dof_satt(fm2_reml)), unname(lmerTest_ddf_fm2), tolerance = 1e-2)
})

## emmeans

if (requireNamespace("emmeans")) {
    salamander1 <- up2date(readRDS(system.file("example_files", "salamander1.rds",
                                                package = "glmmTMB")))

    test_that("emmeans works with ddf='satterthwaite' (GH #1304)", {
        emm <- expect_no_error(
            suppressWarnings(emmeans::emmeans(salamander1, ~ mined, ddf = "satterthwaite"))
        )
        expect_true(all(is.finite(summary(emm)$df)))
    })

    test_that("satterthwaite dffun survives environment-stripping and per-contrast vector calls (GH #1304)", {
        rg <- suppressWarnings(emmeans::ref_grid(salamander1, ddf = "satterthwaite"))

        ## emmeans::ref_grid() replaces dffun's enclosing environment with
        ## baseenv() (glmmTMB#1304); rather than relying on that emmeans
        ## implementation detail, strip it here ourselves so the test keeps
        ## working even if emmeans stops doing so. dffun must not rely on
        ## any free variables (e.g. a captured copy of dof_satt) -- everything
        ## it needs must be reachable via the 'dfargs' argument instead of
        ## lexical scoping
        environment(rg@dffun) <- baseenv()

        ## dffun is called once per contrast with a bare vector (not a
        ## contrast matrix), so the wrapper must turn 'k' into a 1-row
        ## matrix before passing it on to dof_satt()
        k <- rg@linfct[1, ]
        expect_false(is.matrix(k))
        df_val <- suppressMessages(rg@dffun(k, rg@dfargs))
        expect_true(is.finite(df_val))
        expect_length(df_val, 1L)
    })

    ## constructed explicitly via update() (rather than relying on fm1's
    ## REML status, which is mutated earlier in this file if pbkrtest/lme4
    ## are available); update() reuses fm1's existing formula/data, so this
    ## doesn't need direct access to sleepstudy or the lme4 namespace
    fm1_ml_explicit <- update(fm1, REML = FALSE)
    fm1_reml_explicit <- update(fm1_ml_explicit, REML = TRUE)

    test_that("emmeans errors (rather than silently downgrading) for ddf='kenward-roger' on an ML fit", {
        expect_error(
            emmeans::emmeans(fm1_ml_explicit, "Days", ddf = "kenward-roger"),
            "requires a REML fit"
        )
    })

    test_that("emmeans allows ddf='kenward-roger' for a REML fit, with no warning", {
        emm <- expect_no_warning(emmeans::emmeans(fm1_reml_explicit, "Days", ddf = "kenward-roger"))
        expect_true(is.finite(summary(emm)$df[1]))
    })

    test_that("emmeans warns (does not error) for K-R/Satterthwaite on a non-Gaussian family", {
        expect_warning(
            emmeans::emmeans(salamander1, ~ mined, ddf = "satterthwaite"),
            "poorly understood"
        )
    })

    ## salamander1's family (poisson) has no estimated dispersion parameter
    ## (usesDispersion() == FALSE); check_ddf() used to hard-error for this
    ## unconditionally (contradicting the "warn, don't error" GLMM policy
    ## documented above and enforced by emmeans), so summary()/anova()/Anova()
    ## and emmeans() disagreed for exactly this case. satterthwaite's
    ## implementation doesn't need a dispersion parameter and should now warn
    ## (not error) consistently everywhere; kenward-roger's variance-component
    ## machinery genuinely doesn't support such families and should still
    ## error, but with a specific message rather than a crash
    test_that("summary() and emmeans() agree: ddf='satterthwaite' warns (does not error) for a family with no dispersion parameter", {
        expect_warning(summary(salamander1, ddf = "satterthwaite"), "poorly understood")
        expect_warning(emmeans::emmeans(salamander1, ~ mined, ddf = "satterthwaite"), "poorly understood")
    })

    salamander1_reml <- update(salamander1, REML = TRUE)

    test_that("summary() and emmeans() agree: ddf='kenward-roger' errors clearly (not an opaque eigen()/forceSymmetric crash) for a family with no dispersion parameter", {
        expect_error(summary(salamander1_reml, ddf = "kenward-roger"), "no estimated dispersion parameter")
        expect_error(emmeans::emmeans(salamander1_reml, ~ mined, ddf = "kenward-roger"), "no estimated dispersion parameter")
    })

    ## models *without* random effects: emmeans used to decide the ddf on its
    ## own here, rather than deferring to the shared check_ddf() logic that
    ## summary()/anova()/Anova() use
    ss <- transform(lme4::sleepstudy,
                    wt = 1/(1 + Days),
                    grp = factor(Days %% 2))
    fe_free <- glmmTMB(Reaction ~ grp, data = ss)
    fe_pin <- glmmTMB(Reaction ~ grp, data = ss, weights = wt,
                      start = list(betadisp = 0),
                      map = list(betadisp = factor(NA)))
    fe_het <- glmmTMB(Reaction ~ grp, data = ss, dispformula = ~ grp)

    ## df computed from the design rather than read back from the code under
    ## test: 2 fixed-effect coefficients, plus the dispersion parameter when
    ## one is estimated (df.residual.glmmTMB counts it, so these are one
    ## fewer than lm() would report)
    n_obs <- nrow(ss)
    df_free <- n_obs - 3     # 2 beta + betadisp
    df_pin <- n_obs - 2      # betadisp pinned, so not counted

    ## the multiplier emmeans actually used for the confidence limits. This
    ## checks that the df emmeans *reports* is also the one it *used*, which
    ## the df column alone does not; combined with the design-based df above
    ## it pins the value down independently. NB with df = Inf the limits are
    ## named asymp.LCL/asymp.UCL rather than lower.CL/upper.CL
    ci_mult <- function(emm) {
        s <- summary(emm)
        upr <- s[[grep("(UCL|upper\\.CL)$", names(s))[1]]]
        unique(round((upr - s$emmean)/s$SE, 6))
    }

    test_that("emmeans uses asymptotic df when the dispersion parameter is fixed via 'map'", {
        expect_false(glmmTMB:::estDisp(fe_pin))
        expect_equal(sigma(fe_pin), 1)
        expect_equal(df.residual(fe_pin), df_pin)
        ## no variance parameter is estimated at all, so the residual
        ## variance is known and the Wald statistics are exactly standard
        ## normal
        emm <- emmeans::emmeans(fe_pin, ~ grp)
        expect_true(all(is.infinite(summary(emm)$df)))
        expect_equal(ci_mult(emm), round(qnorm(0.975), 6))
    })

    test_that("'map' pinning a dispersion parameter to a shared level still estimates it", {
        ## one parameter is still estimated, so residual df stay in place
        fe_share <- glmmTMB(Reaction ~ grp, data = ss,
                            map = list(betadisp = factor(1)))
        expect_true(glmmTMB:::estDisp(fe_share))
        expect_equal(summary(emmeans::emmeans(fe_share, ~ grp))$df,
                     rep(df_free, 2))
        ## same for a two-parameter dispersion model collapsed to a single
        ## shared level (as opposed to both elements pinned to NA)
        fe_share2 <- glmmTMB(Reaction ~ grp, data = ss, dispformula = ~ grp,
                             map = list(betadisp = factor(c(1, 1))))
        expect_true(glmmTMB:::estDisp(fe_share2))
        fe_pin2 <- glmmTMB(Reaction ~ grp, data = ss, dispformula = ~ grp,
                           start = list(betadisp = c(0, 0)),
                           map = list(betadisp = factor(c(NA, NA))))
        expect_false(glmmTMB:::estDisp(fe_pin2))
    })

    test_that("emmeans respects an explicit ddf='asymptotic' without random effects", {
        ## default is unchanged: residual df, quietly, with the matching
        ## t multiplier for the confidence limits
        expect_no_message(emm <- emmeans::emmeans(fe_free, ~ grp))
        expect_equal(summary(emm)$df, rep(df_free, 2))
        expect_equal(ci_mult(emm), round(qt(0.975, df_free), 6))
        ## an explicit request was previously ignored without any warning
        emm_as <- emmeans::emmeans(fe_free, ~ grp, ddf = "asymptotic")
        expect_true(all(is.infinite(summary(emm_as)$df)))
        expect_equal(ci_mult(emm_as), round(qnorm(0.975), 6))
    })

    test_that("getOption('glmmTMB.df') supplies a default, not an explicit request", {
        ## summary()/anova()/Anova() don't read the option at all, so a value
        ## found there must not override the no-random-effects default -- that
        ## would make emmeans() disagree with them for anyone who has it set,
        ## including someone who sets it to the documented default value
        op <- options(glmmTMB.df = "asymptotic")
        on.exit(options(op), add = TRUE)
        expect_no_message(emm <- emmeans::emmeans(fe_free, ~ grp))
        expect_equal(summary(emm)$df, rep(df_free, 2))
        expect_equal(ci_mult(emm), round(qt(0.975, df_free), 6))
        ## a standing "kenward-roger" must not produce the no-random-effects
        ## message on every call either
        options(glmmTMB.df = "kenward-roger")
        expect_no_message(emm_kr <- emmeans::emmeans(fe_free, ~ grp))
        expect_equal(summary(emm_kr)$df, rep(df_free, 2))
    })

    test_that("emmeans accepts an explicit ddf='df.residual' quietly", {
        ## "df.residual" is what get_ddf() returns internally, and this method
        ## has always accepted it from the user as well. It is deliberately
        ## not one of the match.arg() choices, since summary()/anova()/Anova()
        ## don't take it. Asking for it must not trigger the "no random
        ## effects" message, which is about the K-R/Satterthwaite downgrade
        expect_no_message(emm <- emmeans::emmeans(fe_free, ~ grp,
                                                 ddf = "df.residual"))
        expect_equal(summary(emm)$df, rep(df_free, 2))
        ## and it is honoured even where the new default would be asymptotic:
        ## an explicit request wins over the "nothing estimated" rule
        expect_no_message(emm_pin <- emmeans::emmeans(fe_pin, ~ grp,
                                                     ddf = "df.residual"))
        expect_equal(summary(emm_pin)$df, rep(df_pin, 2))
        expect_equal(ci_mult(emm_pin), round(qt(0.975, df_pin), 6))
    })

    test_that("emmeans rejects an unrecognized ddf instead of falling back to residual df", {
        ## match on the list of choices rather than on match.arg()'s wording,
        ## which is translated (the choices themselves are interpolated
        ## verbatim, so they are locale-independent)
        expect_error(emmeans::emmeans(fe_free, ~ grp, ddf = "no-such-ddf"),
                     "kenward-roger")
        ## the check happens before the component/family branches, so it also
        ## applies to models with random effects (where the error used to come
        ## from check_ddf() further down, with different wording) and to the
        ## zi/disp components (which used to warn and fall back to asymptotic)
        expect_error(emmeans::emmeans(fm1, "Days", ddf = "no-such-ddf"),
                     "kenward-roger")
        expect_error(emmeans::emmeans(fe_het, ~ grp, component = "disp",
                                      ddf = "no-such-ddf"),
                     "kenward-roger")
    })

    test_that("emmeans accepts a partial ddf match, as summary() always has", {
        ## side effect of the new match.arg(): "satt" errored for a model with
        ## random effects before, while summary(ddf = "satt") has always worked
        emm_short <- emmeans::emmeans(fm1, "Days", ddf = "satt")
        emm_full <- emmeans::emmeans(fm1, "Days", ddf = "satterthwaite")
        expect_equal(summary(emm_short)$df, summary(emm_full)$df)
        expect_false(any(is.infinite(summary(emm_short)$df)))
    })

    test_that("ddf for a component other than 'cond' is still asymptotic", {
        ## unchanged behaviour, but it now passes through the new match.arg()
        expect_true(all(is.infinite(summary(emmeans::emmeans(
            fe_het, ~ grp, component = "disp", ddf = "asymptotic"))$df)))
        expect_warning(emm <- emmeans::emmeans(fe_het, ~ grp,
                                               component = "disp",
                                               ddf = "satterthwaite"),
                       "using ddf 'asymptotic' instead")
        expect_true(all(is.infinite(summary(emm)$df)))
    })

    ## the two blocks below pin down emmeans behaviour that must *not* change:
    ## their emmeans expectations hold on unpatched master too (the estDisp()
    ## assertion below is new, since the helper itself is)
    test_that("emmeans default for a Gaussian model with dispformula = ~0 is asymptotic", {
        fe_zero <- glmmTMB(Reaction ~ grp, data = ss, dispformula = ~ 0)
        expect_false(glmmTMB:::estDisp(fe_zero))
        expect_true(all(is.infinite(summary(emmeans::emmeans(fe_zero, ~ grp))$df)))
    })

    test_that("emmeans behaviour for a non-Gaussian model without random effects is unchanged", {
        ss$cnt <- round(ss$Reaction)
        fe_pois <- glmmTMB(cnt ~ grp, family = poisson, data = ss)
        ## asymptotic by default ...
        expect_true(all(is.infinite(summary(emmeans::emmeans(fe_pois, ~ grp))$df)))
        ## ... and an explicit kenward-roger/satterthwaite request is still
        ## downgraded to asymptotic with a warning, before the no-random-effects
        ## logic above is reached
        for (dd in c("satterthwaite", "kenward-roger")) {
            expect_warning(emm <- emmeans::emmeans(fe_pois, ~ grp, ddf = dd),
                           "using ddf 'asymptotic' instead")
            expect_true(all(is.infinite(summary(emm)$df)))
        }
    })

    test_that("summary() and emmeans() agree: K-R/Satterthwaite fall back to residual df without random effects", {
        ## df from the design: 2 fixed-effect coefficients plus however many
        ## dispersion parameters are estimated
        mods <- list(fe_free, fe_het, fe_pin)
        exp_df <- c(df_free, n_obs - 4, df_pin)
        for (dd in c("satterthwaite", "kenward-roger")) {
            for (i in seq_along(mods)) {
                ## dof_satt()/dof_KR() need random effects; for a non-trivial
                ## dispformula emmeans used to reach them anyway and error,
                ## and for a trivial one it silently ignored the request
                expect_message(emm <- emmeans::emmeans(mods[[i]], ~ grp, ddf = dd),
                               "no random effects in model")
                expect_equal(summary(emm)$df, rep(exp_df[i], 2))
                expect_equal(ci_mult(emm), round(qt(0.975, exp_df[i]), 6))
                ## ... and that is what summary() reports too
                expect_equal(
                    unname(suppressMessages(
                        summary(mods[[i]], ddf = dd)$coefficients$cond[, "ddf"])),
                    rep(exp_df[i], 2))
            }
        }
    })

    test_that("satterthwaite ddf works for a rank-deficient conditional model (#1326)", {
        rd_data <- data.frame(
            genotype = factor(rep(c("a", "c", "b", "c"), c(6L, 6L, 12L, 6L))),
            plate = factor(rep(c("1", "2", "3", "4", "5", "6"), 5)),
            n_seeds = 5,
            n_disease = rep(
                c(0, 1, 3, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0),
                c(1L, 1L, 1L, 6L, 1L, 2L, 2L, 2L, 1L, 7L, 1L, 2L, 1L, 2L)
            ),
            year = factor(rep(c("1", "2"), c(18L, 12L)))
        )
        fit_rd <- expect_message(
            glmmTMB(cbind(n_disease, n_seeds - n_disease) ~ genotype * year +
                        (1 | plate:genotype:year),
                    family = binomial, data = rd_data),
            "dropping columns from rank-deficient conditional model"
        )

        ## the aliased coefficient (genotypec:year2) gets NA ddf, like its
        ## NA estimate/SE, rather than erroring
        cc <- suppressWarnings(summary(fit_rd, ddf = "satterthwaite"))$coefficients$cond
        expect_true(is.na(cc["genotypec:year2", "ddf"]))
        expect_true(all(is.finite(cc[rownames(cc) != "genotypec:year2", "ddf"])))

        emm <- suppressWarnings(emmeans::emmeans(
            fit_rd, ~ genotype | year, weights = "flat",
            nesting = "genotype %in% year", ddf = "satterthwaite"))
        expect_true(all(is.finite(summary(emm)$df)))
    })
}

## ---- coefficients tied/fixed via 'map', with and without rank deficiency (GH #1340) ----
##
## Oracles: a coefficient fixed to zero via 'map' is the same model as the
## formula without that column, and coefficients tied via 'map' are the same
## model as the formula with the columns summed; both refits go through the
## unmapped code path. The ddf of the mapped fit must reproduce them.
if (requireNamespace("lme4")) {
    ss <- lme4::sleepstudy
    ss$Days2 <- ss$Days^2
    f_tied <- glmmTMB(Reaction ~ Days + Days2 + (1 | Subject), data = ss, REML = TRUE,
                      map = list(beta = factor(c(1, 2, 2))))
    f_fixed <- glmmTMB(Reaction ~ Days + Days2 + (1 | Subject), data = ss, REML = TRUE,
                       start = list(beta = c(250, 10, 0)),
                       map = list(beta = factor(c(1, 2, NA))))
    o_tied <- glmmTMB(Reaction ~ I(Days + Days2) + (1 | Subject), data = ss, REML = TRUE)
    o_fixed <- glmmTMB(Reaction ~ Days + (1 | Subject), data = ss, REML = TRUE)
    ## (the Satterthwaite tolerance allows for its finite-difference Jacobians)
    tol_kr <- 1e-8
    tol_satt <- 1e-6

    test_that("mapped fits are the same models as their unmapped oracles", {
        expect_equal(logLik(f_tied), logLik(o_tied), tolerance = 1e-8, ignore_attr = TRUE)
        expect_equal(logLik(f_fixed), logLik(o_fixed), tolerance = 1e-8, ignore_attr = TRUE)
    })

    test_that("ddf for coefficients tied via 'map' (one value per coefficient, the tied ones shared)", {
        for (ff in list(list(dof_KR, tol_kr), list(dof_satt, tol_satt))) {
            d <- unname(c(ff[[1]](f_tied)))
            expect_length(d, 3L)
            expect_identical(d[2], d[3])
            expect_equal(d[1:2], unname(c(ff[[1]](o_tied))), tolerance = ff[[2]])
        }
        ## the KR-adjusted vcov is reported for all three coefficients, with
        ## the tied pair perfectly correlated
        V <- attr(dof_KR(f_tied), "vcov")
        expect_identical(dim(V), c(3L, 3L))
        expect_equal(V["Days", ], V["Days2", ], ignore_attr = TRUE)
        expect_equal(unname(V[1:2, 1:2]), unname(attr(dof_KR(o_tied), "vcov")), tolerance = tol_kr)
    })

    test_that("ddf for a coefficient fixed via 'map' is NA (known constant), the others unchanged", {
        for (ff in list(list(dof_KR, tol_kr), list(dof_satt, tol_satt))) {
            d <- unname(c(ff[[1]](f_fixed)))
            expect_length(d, 3L)
            expect_true(is.na(d[3]))
            expect_equal(d[1:2], unname(c(ff[[1]](o_fixed))), tolerance = ff[[2]])
        }
        ## zero (not NA) variance in the reported KR vcov, so downstream
        ## linear algebra can use it; NA in the printed summary
        V <- attr(dof_KR(f_fixed), "vcov")
        expect_identical(dim(V), c(3L, 3L))
        expect_equal(unname(V["Days2", ]), rep(0, 3))
        cc <- summary(f_fixed, ddf = "kenward-roger")$coefficients$cond
        expect_true(is.na(cc["Days2", "ddf"]))
        expect_identical(unname(cc[, "ddf"]), unname(c(dof_KR(f_fixed))))
    })

    test_that("Satterthwaite ddf of a contrast: NA only if it involves *only* map-fixed coefficients", {
        L <- rbind(c(0, 0, 1),   # Days2 alone: fixed, no variance
                   c(0, 1, 1))   # Days + Days2: estimable, same as Days
        d <- dof_satt(f_fixed, L = L)
        expect_true(is.na(d[1]))
        expect_equal(d[2], unname(dof_satt(o_fixed))[2], tolerance = tol_satt)
    })

    ## empty cell c:2 -> the g:yr interaction loses a column
    ss2 <- transform(ss,
                     g = factor(rep(c("a", "b", "c"), length.out = nrow(ss))),
                     yr = factor(ifelse(ss$Days < 5, "1", "2")))
    ss2$g[ss2$g == "c" & ss2$yr == "2"] <- "b"
    ss2$g <- droplevels(ss2$g)

    test_that("ddf for a rank-deficient model with tied coefficients (both reductions at once)", {
        ## on top of the dropped column, the two g main-effect coefficients are tied
        f_both <- suppressMessages(
            glmmTMB(Reaction ~ g * yr + (1 | Subject), data = ss2, REML = TRUE,
                    map = list(beta = factor(c(1, 2, 2, 3, 4)))))
        expect_length(fixef(f_both)$cond, 6L)          # nominal
        expect_identical(ncol(getME(f_both, "X")), 5L)   # rank-reduced
        ## oracle: refit on the estimated parameters' design columns
        X <- as.matrix(getME(f_both, "X"))
        Z <- cbind(X[, 1], X[, 2] + X[, 3], X[, 4], X[, 5])
        colnames(Z) <- paste0("Z", 1:4)
        o_both <- glmmTMB(Reaction ~ 0 + Z1 + Z2 + Z3 + Z4 + (1 | Subject),
                          data = data.frame(Reaction = ss2$Reaction, Subject = ss2$Subject, Z),
                          REML = TRUE)
        expect_equal(logLik(f_both), logLik(o_both), tolerance = 1e-8, ignore_attr = TRUE)
        for (ff in list(list(dof_KR, tol_kr), list(dof_satt, tol_satt))) {
            d <- unname(c(ff[[1]](f_both)))
            o <- unname(c(ff[[1]](o_both)))
            expect_length(d, 6L)
            expect_true(is.na(d[6]))                    # aliased
            expect_equal(d[c(1, 2, 4, 5)], o, tolerance = ff[[2]])
            expect_identical(d[2], d[3])                # tied
        }
    })

    test_that("a map that only permutes the parameters gives the ddf of the unmapped fit", {
        f_perm <- glmmTMB(Reaction ~ Days + Days2 + (1 | Subject), data = ss, REML = TRUE,
                          map = list(beta = factor(c(2, 1, 3))))
        f_plain <- glmmTMB(Reaction ~ Days + Days2 + (1 | Subject), data = ss, REML = TRUE)
        expect_equal(dof_KR(f_perm), dof_KR(f_plain), tolerance = tol_kr, ignore_attr = TRUE)
        expect_equal(dof_satt(f_perm), dof_satt(f_plain), tolerance = tol_satt)
    })

    test_that("anova() F tests carry the restriction over to the tied parameter", {
        m0 <- glmmTMB(Reaction ~ 1 + (1 | Subject), data = ss, REML = TRUE)
        for (dd in c("kenward-roger", "satterthwaite")) {
            a <- anova(m0, f_tied, ddf = dd)
            o <- anova(m0, o_tied, ddf = dd)
            expect_equal(a[2, "Num Df"], 1)      # two tied columns, one restriction
            expect_equal(unlist(a[2, c("F", "Num Df", "Den Df")]),
                         unlist(o[2, c("F", "Num Df", "Den Df")]), tolerance = tol_satt)
        }
    })

    test_that("a joint hypothesis that contradicts a nonzero fixed value gives NA with a warning", {
        f_half <- glmmTMB(Reaction ~ Days + Days2 + (1 | Subject), data = ss, REML = TRUE,
                          start = list(beta = c(250, 10, 0.5)), map = list(beta = factor(c(1, 2, NA))))
        L <- rbind(c(0, 1, 0), c(0, 1, 1))  # Days = 0 and Days + Days2 = 0, but Days2 == 0.5
        expect_warning(res <- glmmTMB:::.joint_test(f_half, L, "kenward-roger"), "contradict")
        expect_true(is.na(res$Fstat))
        ## with the coefficient fixed to zero the two rows agree and reduce to one
        res0 <- glmmTMB:::.joint_test(f_fixed, L, "kenward-roger")
        expect_equal(res0$ndf, 1)
        expect_equal(res0$Fstat, glmmTMB:::.joint_test(f_fixed, L[1, , drop = FALSE], "kenward-roger")$Fstat)
    })

    test_that("ddf are NA (without error) when every coefficient is fixed via 'map'", {
        f_all <- glmmTMB(Reaction ~ 1 + (1 | Subject), data = ss, REML = TRUE,
                         start = list(beta = 298), map = list(beta = factor(NA)))
        expect_true(is.na(dof_KR(f_all)))
        expect_true(is.na(dof_satt(f_all)))
    })

    test_that("an NA-filled fixed-effect vcov (boundary fit) gives NA ddf with a warning, not an eigen() error", {
        ## sdreport() can return one for a fit on the boundary; inject it
        testthat::with_mocked_bindings(
            .Phi_est = function(model, component = "cond") matrix(NA_real_, 2, 2),
            .package = "glmmTMB", {
                expect_warning(d <- dof_KR(o_fixed), "contains NA values")
                expect_true(all(is.na(d)))
                expect_warning(d <- dof_satt(o_fixed), "contains NA values")
                expect_true(all(is.na(d)))
            })
    })

    if (requireNamespace("emmeans")) {
        test_that("emmeans works for coefficients tied via 'map' (all ddf), matching the unmapped oracle", {
            for (dd in c("asymptotic", "kenward-roger", "satterthwaite")) {
                em <- summary(emmeans::emmeans(f_tied, ~ 1, at = list(Days = 3, Days2 = 9), ddf = dd))
                eo <- summary(emmeans::emmeans(o_tied, ~ 1, at = list(Days = 3, Days2 = 9), ddf = dd))
                expect_equal(em$emmean, eo$emmean, tolerance = 1e-8)
                expect_equal(em$SE, eo$SE, tolerance = 1e-8)
                expect_equal(em$df, eo$df, tolerance = tol_satt)
            }
        })
        test_that("emmeans works for a coefficient fixed via 'map' with ddf='kenward-roger'", {
            em <- summary(emmeans::emmeans(f_fixed, ~ 1, at = list(Days = 3, Days2 = 9), ddf = "kenward-roger"))
            eo <- summary(emmeans::emmeans(o_fixed, ~ 1, at = list(Days = 3), ddf = "kenward-roger"))
            expect_equal(em$SE, eo$SE, tolerance = 1e-8)
            expect_equal(em$df, eo$df, tolerance = tol_kr)
        })
    }

    if (requireNamespace("car")) {
        test_that("car::Anova F tests run for a rank-deficient model", {
            f_rd <- suppressMessages(glmmTMB(Reaction ~ g * yr + (1 | Subject), data = ss2, REML = TRUE))
            w <- car::Anova(f_rd)   # Type II Wald, on the non-aliased coefficients
            for (dd in c("kenward-roger", "satterthwaite")) {
                a <- car::Anova(f_rd, ddf = dd)
                expect_equal(a[["Num Df"]], w[["Df"]])
                expect_true(all(is.finite(a[["F"]])))
            }
            ## Satterthwaite only changes the denominator: F * Num Df is the Wald chi-square
            expect_equal(a[["F"]] * a[["Num Df"]], w[["Chisq"]], tolerance = 1e-8)
        })
        test_that("car::Anova F tests work for coefficients tied or fixed via 'map'", {
            for (dd in c("kenward-roger", "satterthwaite")) {
                ## tied: each of Days/Days2 tests the one shared parameter
                a <- car::Anova(f_tied, type = "III", ddf = dd)
                o <- car::Anova(o_tied, type = "III", ddf = dd)
                expect_equal(unlist(a["Days", ]), unlist(a["Days2", ]))
                expect_equal(unlist(a["Days", ]), unlist(o["I(Days + Days2)", ]), tolerance = tol_satt)
                ## fixed: an untestable (NA) row for the fixed coefficient
                a <- car::Anova(f_fixed, type = "III", ddf = dd)
                o <- car::Anova(o_fixed, type = "III", ddf = dd)
                expect_true(is.na(a["Days2", "F"]))
                expect_equal(unlist(a["Days", ]), unlist(o["Days", ]), tolerance = tol_satt)
            }
        })
    }
}
