# Leverage / hatvalues (PR #1065): status and to-do

## Where things stand (2026-11-23)

- `hatvalues(model, type = "full")` is Kasper's three-path calculation. Checked against brute-force refits for rr nbinom2, Poisson GLMM and LMM.
- `hatvalues(model, type = "conditional")` holds covariance/dispersion parameters fixed and treats `beta`/`betazi` jointly with the random effects (the `lme4` definition). It needs no `newton()` and no third derivatives.
  - Matches `lme4::hatvalues()` for LMMs to about 1e-6.
  - For GLMMs it is within about 1e-4 of the exact "re-estimate beta only, via the Laplace marginal likelihood" leverage.
  - Memory for an rr nbinom2 model with nobs = 8000: 556 MB and 1.7 s.

## Diagnosis of the memory blow-up for rr models

- All of the cost is in `F$laplace(r)$newton(p)` (path 1, y -> theta_hat). `F$laplace(r)` alone costs almost nothing.
- Memory grows as about 68 KB x nobs x ntheta for rr nbinom2, and about 30 KB x nobs x ntheta for a `diag()` control with the same number of random effects.
  - For rr, ntheta is about ntax*d, so memory grows roughly quadratically in the number of taxa: about 200 GB for ntax = 200, nsubj = 40, d = 2.
- The rr Hessian of the random effects (Hbb) is not dense: it is block-diagonal with d x d blocks per subject. The non-spherical rr reparameterization therefore wouldn't help.
- `TMB::config(tmbad.atomic_sparse_log_determinant = 1)` fails in `newton()` with "Inverse subset: order 2 not yet implemented".
- Retaping with `ADreport = TRUE` also builds `ADGrad` for the nobs-output prediction function unless `env$type` is restricted to `"ADFun"`. That costs several GB on large models. The conditional code restricts it; the full code does not yet.

## To do

### Testing
- [ ] Turn the ad hoc checks into `tests/testthat` tests, skipped if RTMB isn't installed:
  - conditional vs `lme4::hatvalues` (LMM);
  - no-RE model vs `lm`;
  - `diag = FALSE` diagonal vs `diag = TRUE`;
  - a small brute-force comparison (full and conditional).
- [ ] Brute-force step size: use eps of about 1e-2. With eps = 1e-4 the refit stops at its starting optimum and gives wrong leverages, e.g. identical values for different observations in an LMM. `leverage_brute()` (misc/leverage_funs.R) defaults to 1e-3, so check earlier comparisons made with it.
- [ ] More families and model types:
  - binomial with n > 1, Beta, Gamma, Tweedie, truncated families;
  - `dispformula` with random effects;
  - zero-inflation with random effects (works; negative hat values for y = 0 appear in both versions and look genuine, but need a sanity check);
  - `weights`, offsets, NA handling (`na.exclude`), `map`ped parameters.
- [ ] Discrete responses: decide whether d mu / d y is the right definition. The derivative is taken as though y were continuous.
- [ ] Scaling of the conditional version beyond nobs = 8000 (e.g. ntax = 400). Also models with many fixed effects (e.g. `0 + taxon` intercepts), where the `beta` row of Hxx is dense.

### Code
- [ ] Restrict `env$type` to `"ADFun"` in the full version too (path 2 retapes with `ADreport = TRUE`).
- [ ] Full version for rr models: restrict path 1's `newton()` to a subset of parameters, or compute J by the implicit function theorem from the fitted Hessian of the fixed effects.
- [ ] Remove the RTMB dependency (`GetTape`, `MakeTape`, `advector`, `TMB::config`). The conditional version only needs a sparse Hessian over (yobs, x) and a sparse Jacobian of mu, so it may be the easiest to port to plain TMB.
- [ ] `TMB:::solveSubset` is private: export it from TMB or copy it.
- [ ] Compiling without OpenMP is needed for RTMB compatibility, and `-DTMB_MAX_ORDER=4` is needed for `newton()`. Both changes to `src/Makevars` must be resolved before merging.
- [ ] rr terms allocate ntax x nsubj spherical u's, but only d x nsubj enter the linear predictor. Map the rest off, or allocate rank x blockReps (separate PR: touches startVals, ranef, simulate, fix_simcode).
- [ ] `rstudent.glmmTMB`: account for sigma(-i), handle zero weights, `na.exclude` and naming.
- [ ] Decide the default `type`. "full" is the more complete definition; "conditional" is much cheaper and matches `lme4`.
