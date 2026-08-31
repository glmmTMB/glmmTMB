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
}
