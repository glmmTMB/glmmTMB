# Satterthwaite ddf performance: findings and changes (`satterthwaite_modifications` branch)

## Background

While running a Type-I-error simulation comparing `lmerTest`, `nlme::lme`,
and `glmmTMB` (using the `ddf`-corrected `anova()`/`summary()` from PR #1302,
branch `anova_ddf`) on heteroscedastic data, the Satterthwaite/Kenward-Roger
batch of the simulation was ~10x slower per replicate than every other
method combined. This document records what turned out to be causing that,
what was fixed here (on `master`, independent of the still-unmerged PR), and
what's left for later.

This branch is **based on `master`**, not on `anova_ddf`, so it only touches
`dof_satt()` (used by `summary.glmmTMB(object, ddf = "satterthwaite")`),
which already existed on master. The PR's `.joint_ddf_satt()` (used by
`anova.glmmTMB(..., ddf = "satterthwaite")`) doesn't exist here — see
"Future work" below for how to extend this to it.

## Diagnosis

Benchmarking a single replicate of the simulation's `glmmTMB` REML fit
(`y ~ trt + (1|id)`, `dispformula = ~trt`) broke down like this:

| component | mean time | % of total |
|---|---|---|
| fit full model (REML) | 0.063s | 5.0% |
| fit null model (REML) | 0.058s | 4.6% |
| `anova(ddf="satterthwaite")` | 0.515s | 40.5% |
| `anova(ddf="kenward-roger")` | 0.057s | 4.5% |
| `summary(ddf="satterthwaite")` | 0.514s | 40.4% |
| `summary(ddf="kenward-roger")` | 0.065s | 5.1% |

Two things stood out:

1. **Model fitting itself is cheap (~9.5%); the ddf calculation dominates
   (~90.5%).**
2. **Satterthwaite is ~9x more expensive than Kenward-Roger**, despite KR
   usually having the reputation for being the expensive correction.

### Why Satterthwaite is so much slower than Kenward-Roger

`dof_KR()`'s path (`.vcov_kenward_adjusted()` → `.vcovAdj16_internal()` /
`.get_SigmaG()`, in `denom_df.R`) is pure closed-form linear algebra: it
builds the classical pbkrtest P/W sensitivity matrices analytically from
`vcov(model)$cond` and the design matrix. No numerical differentiation, no
re-evaluation of the model at perturbed parameter values.

`dof_satt()`'s path, by contrast, computes the same kind of sensitivity via
**nested finite differencing**:

```r
kappa_opt <- model$fit$par
h_kappa   <- numDeriv::jacobian(func = model$obj$gr, x = kappa_opt)   # ~2p evaluations of obj$gr
jac_kappa <- .get_jac_list(.covbeta_kappa, kappa_opt, model)          # ~2p evaluations of .covbeta_kappa
```

and `.covbeta_kappa()` does, *per perturbed point*:

```r
sdr <- TMB::sdreport(md$obj, par.fixed = kappa, getJointPrecision = TRUE)
```

`sdreport()` is a heavyweight, general-purpose routine (assembles the full
joint precision matrix, does a GMRF marginalization, computes delta-method
quantities for every reported value) — and it was being called roughly `2p`
times (`p` = number of variance/dispersion parameters) just to numerically
differentiate one covariance matrix.

**Caveat on the KR/Satterthwaite comparison**: it isn't entirely apples to
apples. `.get_SigmaG()` folds the residual variance in as a single scalar —
the classical, homoscedastic pbkrtest assumption — which matches the message
glmmTMB prints for non-trivial `dispformula` models:
`"ddf='kenward-roger' ignored except for conditional-distribution
parameters"`. KR's ddf calculation does **not** propagate dispersion-formula
parameter uncertainty into the correction. Satterthwaite's `kappa <-
model$fit$par` includes *all* parameters (beta, theta, and the dispersion
regression coefficients), so it's correctly accounting for more sources of
uncertainty for exactly the heteroscedastic-dispersion models this
simulation cares about — which is also part of why it does more work. Any
future attempt to "port KR's cheap approach to Satterthwaite" needs to keep
this in mind, or it would silently drop that correctness for models with
non-trivial `dispformula`.

### Comparison with lmerTest's Satterthwaite implementation

`lmerTest` (`~/R/pkgs/lmerTestR`) solves the same statistical problem far
more cheaply, for two independent reasons:

**1. It caches the expensive pieces once, at fit time.** In
`as_lmerModLT()` (`R/lmer.R`), called once when you fit via
`lmerTest::lmer()`:

```r
h <- numDeriv::hessian(func = devfun_vp, x = varpar_opt, devfun = devfun, reml = is_reml)
...
res@vcov_varpar <- 2 * h_inv
Jac <- numDeriv::jacobian(func = get_covbeta, x = varpar_opt, devfun = devfun)
res@Jac_list <- lapply(...)
```

Every subsequent `summary()`/`anova()`/`contest()` call just reads
`@vcov_varpar`/`@Jac_list` off the object. Zero recomputation, ever, for a
given fit.

**2. Its inner function is cheap.** lmerTest's Jacobian is of `get_covbeta()`:

```r
get_covbeta <- function(varpar, devfun) {
  sigma <- varpar[nvarpar]; theta <- varpar[-nvarpar]
  devfun(theta)                                        # one profiled PLS/Cholesky update
  sigma^2 * tcrossprod(environment(devfun)$pp$RXi())   # closed-form GLS cov(beta)
}
```

This is cheap because lme4's `devfun` is purpose-built to let you cheaply
re-profile beta given theta via a single sparse-Cholesky update — no
re-optimization, no Laplace approximation (a Gaussian LMM's marginal
likelihood is exact, so there's no Laplace error to approximate away in the
first place). glmmTMB's `.covbeta_kappa()` instead reaches for the fully
general `TMB::sdreport()`, built to handle arbitrary non-Gaussian,
Laplace-approximated models, which does far more work per evaluation than
lmerTest's one-line GLS extraction.

## What was changed here

Two changes to `dof_satt()` in `glmmTMB/R/denom_df.R`, both **implemented,
no algorithmic change to the Satterthwaite formula itself**:

### 1. Cache the per-model pieces (`.satt_precompute()`)

`kappa_opt`, `cov_varpar_kappa`, and `jac_kappa` depend only on the fitted
model, not on the contrast `L` being tested. They're now computed once by a
new helper, `.satt_precompute(model)`, and cached on `model$obj$env` — a
genuine R environment (unlike the rest of `model`, an ordinary list), so it
is shared by reference across every copy of `model` and persists across
separate top-level calls (e.g. two separate `summary(fit, ddf =
"satterthwaite")` calls on the same fit):

```r
.satt_precompute <- function(model) {
    cache_env <- model$obj$env
    if (!is.null(cache_env$.satt_cache)) {
        return(cache_env$.satt_cache)
    }
    kappa_opt <- model$fit$par
    h_kappa <- numDeriv::jacobian(func = model$obj$gr, x = kappa_opt, method = "simple")
    eig_h_kappa <- eigen(h_kappa, symmetric = TRUE)
    cov_varpar_kappa <- with(eig_h_kappa, vectors %*% diag(1 / values) %*% t(vectors))
    jac_kappa <- .get_jac_list(.covbeta_kappa, kappa_opt, model, method = "simple")
    res <- list(cov_varpar_kappa = cov_varpar_kappa, jac_kappa = jac_kappa)
    cache_env$.satt_cache <- res
    res
}
```

`dof_satt()` now calls this instead of inlining the computation.

### 2. Cheaper finite differencing (`method = "simple"`)

Both `numDeriv::jacobian()` calls previously used the default
`method = "Richardson"`, which re-evaluates each finite difference at
several step sizes for extrapolated accuracy (~4-8x the function
evaluations of a one-sided difference). Since the end result is an
approximate denominator df (not a precision-critical quantity), both calls
now use `method = "simple"` (~p+1 evaluations total). A prior, already
-present comment in the code (`## does not help anything -- but seems
precision is already good enough: method.args = list(r = 6)`) recorded that
someone had already checked whether *more* Richardson precision mattered and
found it didn't — consistent with there being little to lose by using
*less*-precise (and much cheaper) differencing here.

## Results

Benchmarked on 15 replicates of the simulation's exact `glmmTMB` REML fit
(`y ~ trt + (1|id)`, `dispformula = ~trt`), comparing the old inline
Richardson/no-cache code against the new cached/simple-differencing
`dof_satt()`:

```
---- ddf agreement (old vs new, method=Richardson vs simple) ----
absolute difference:  mean 0.0082  (range 0.0051 - 0.0132; ddf values are ~20-40)
relative difference:  mean 0.026%  (range 0.023% - 0.027%)

---- timing ----
old (Richardson, no cache): mean 0.475s/call
new (simple, cached):       mean 0.082s/call
speedup: 5.8x

repeated dof_satt() on same model: 1st call 0.080s, 2nd call 0.0129s (cached)
```

- Denominator df values agree to within ~0.03% — negligible for an
  already-approximate quantity, and p-values from
  `summary(fit, ddf="satterthwaite")` matched to 4+ decimal places
  end-to-end.
- A single (uncached) call is **5.8x faster** from the differencing-method
  change alone.
- A **second** call on the same fitted model (the caching win) drops to
  ~0.013s — roughly another 6x on top, ~36x faster than the original for
  repeat calls.
- `dof_KR()` / `ddf = "kenward-roger"` is untouched and gives identical
  results to before, as expected.

## Test coverage added

`tests/testthat/test-ddf.R` previously had no test pinning `dof_satt()`'s
Satterthwaite values against an independent external reference (unlike
`dof_KR()`, which is checked against `pbkrtest` directly) — only a
self-consistency check that `summary(fit, ddf="satterthwaite")` matches
`dof_satt(fit)`, which would pass even if both were wrong in the same way.
Added `test_that("Satterthwaite ddf match lmerTest (hard-coded reference
values)", ...)`, comparing `dof_satt()` on the existing `fm1`/`fm2` sleepstudy
fixtures against `lmerTest::lmer()`'s Satterthwaite df for the same models
(REML, same formulas). The lmerTest values are hard-coded in the test rather
than computed at test time, to avoid adding a test dependency on lmerTest.
Tolerance is 1%, comfortably above the ~0.03% difference actually observed
between glmmTMB's (simple-differencing) and lmerTest's implementations,
while still tight enough to catch a real regression.

## Future work

Not implemented here; recorded for later:

1. **Wire this into the `anova_ddf` PR branch.** The redundancy that
   originally motivated this (a single simulation replicate calling both
   `anova(null, full, ddf="satterthwaite")` *and*
   `summary(full, ddf="satterthwaite")`, each independently paying the full
   cost) involves `.joint_ddf_satt()`, which only exists on `anova_ddf`, not
   here. Once this branch and `anova_ddf` are combined, `.joint_ddf_satt()`
   should be refactored to call `.satt_precompute()` too, so the two code
   paths share the cache. That's the change that will actually eliminate the
   specific redundancy the simulation hit (summary + anova on the same
   model), rather than just speeding up repeated calls to `dof_satt()` alone.

2. **Replace `.covbeta_kappa()`'s `TMB::sdreport()` call with something
   narrower.** It only needs one marginal block (`Cov(beta | kappa)`) of the
   joint precision, but `sdreport(getJointPrecision = TRUE)` computes and
   packages far more than that. Calling TMB's lower-level internals directly
   to get just the joint Hessian and doing the `GMRFmarginal()` step by hand
   should speed up each of the ~p+1 remaining calls, on top of the fixes
   here.

3. **Gaussian-family-specific closed-form shortcut.** For the Gaussian
   family specifically (which covers this simulation), beta given the rest
   of the parameters is *still* available in closed form even with
   heteroscedastic dispersion — it's a generalized least squares solve with
   `V = Z G Z' + R(kappa)`, no Laplace approximation needed (Gaussian
   marginal likelihoods are exact, unlike the general Laplace-approximated
   case `sdreport()` is built for). A `get_covbeta_glmmTMB(kappa)` analogous
   to lmerTest's `get_covbeta()` — assembling `V` directly and solving the
   GLS normal equations — could replace `sdreport()` for Gaussian-family
   models entirely, closing most of the remaining gap to lmerTest's speed
   while (unlike Kenward-Roger's shortcut) still correctly including
   dispersion-parameter uncertainty. Would not generalize to non-Gaussian
   families (Poisson/binomial/etc.), where lme4 has no shortcut either and
   the full Laplace/`sdreport()` machinery is genuinely necessary.

4. **Analytic (AD-based) sensitivity instead of any finite differencing.**
   Since TMB already computes `Cov(beta)` via automatic differentiation, its
   derivative with respect to `kappa` could in principle be obtained via AD
   (implicit-function-theorem style) rather than finite-differencing an
   expensive function at all. This is the "correct" long-term fix but a much
   bigger implementation lift than 1-3, and worth its own discussion with
   the TMB maintainers before attempting.

5. **Parallelize the finite-difference evaluations**, for the (different)
   use case of a single one-off large model rather than many parallel
   simulation replicates: `numDeriv::jacobian`'s evaluations at different
   perturbed points are embarrassingly parallel and could be farmed out with
   `future.apply`/`parallel::mclapply` if there are idle cores within a
   single fit. Not useful for a many-replicates simulation where all cores
   are already busy across replicates (as in the simulation that prompted
   this investigation).

## Files changed

- `glmmTMB/R/denom_df.R`: added `.satt_precompute()`; `dof_satt()` now uses
  it instead of inlining the Hessian/Jacobian computation.
- `glmmTMB/tests/testthat/test-ddf.R`: added a Satterthwaite-vs-lmerTest
  regression test (hard-coded reference values, no new test dependency).

No other files were touched; `dof_KR()`/Kenward-Roger and the C++ side are
unaffected. Existing tests: all pass (`test-ddf.R`: 6/6 after the addition
above; no other test file in the suite touches `ddf`/Satterthwaite/KR).
