# Separable/Kronecker Prototype Plan

This is a local working plan for replacing the specialized `xar1` implementation
with a more generic separable covariance backend. It should stay untracked unless
it is later converted into package documentation or developer notes.

## Goal

Prototype a generic internal separable random-effect covariance engine that can
reproduce the two current structures:

```r
homcsxar1(...)
unxar1(...)
```

as separable products:

```r
separable(sepgrid(member, time) + 0 | group,
          homcs(member),
          ar1(time))

separable(sepgrid(member, time) + 0 | group,
          us(member),
          ar1(time))
```

The immediate goal is not a fully general "any structure x any structure"
interface. Additional combinations should be added deliberately because margin
structures differ in parameter counts, coordinate requirements, scale
conventions, and efficient TMB representations.

## Scope For First Prototype

- Add one internal `kron_covstruct` or `separable_covstruct` code path.
- Support two margin combinations:
  - `homcs(member) x ar1(time)`
  - `us(member) x ar1(time)`
- Keep scale in the member margin.
- Treat the AR(1) time margin as correlation-only.
- Require a full rectangular latent grid.
- Allow missing observed rows when factor/coordinate levels define the intended
  full grid.
- Keep `homcsxar1()` and `unxar1()` as aliases or wrappers if maintainers want
  the dyad-friendly names.

## Syntax Strategy

Start with low-level explicit syntax:

```r
separable(sepgrid(member, time) + 0 | group,
          homcs(member),
          ar1(time))
```

This keeps `| group` in the usual random-effect position and makes the latent
grid explicit. Later, after the backend works, consider higher-level sugar such
as:

```r
kron(homcs(member), ar1(time), group = group)
```

or:

```r
separable(homcs(member), ar1(time), group = group)
```

## R-Side Tasks

1. Add parser support for `separable(<grid + 0 | group>, <margin1>, <margin2>)`.
2. Add `sepgrid()` as a general Cartesian coordinate helper.
3. Validate that:
   - the grid has exactly two coordinates for the first prototype;
   - grid variables match margin variables;
   - `+ 0` is used;
   - levels define a complete Cartesian product;
   - the AR(1) time coordinate is discrete and unit-spaced;
   - repeated coordinate variables are rejected unless explicitly supported.
4. Compute and store per-term metadata:
   - separable dimensions;
   - margin structure codes;
   - theta lengths;
   - scale-carrying margin;
   - coordinate levels;
   - permutation between model-matrix column order and array order.
5. Keep output labels stable and compact.

## C++ / TMB Tasks

1. Add a separable branch in `termwise_nll()`.
2. Split theta by margin.
3. Build member-margin covariance objects:
   - homogeneous compound symmetry;
   - unstructured covariance.
4. Transform the AR(1) parameter using the existing `glmmTMB` convention.
5. Reshape each group-specific random-effect block into a `member x time` array.
6. Evaluate the separable prior without constructing the full Kronecker matrix.
   Benchmark both:
   - `SEPARABLE(AR1(phi), member_density)(x)`;
   - `AR1(phi, member_density)(x)`.
7. Manually apply member-specific scale and log-Jacobian if needed.
8. Add simulation support only after likelihood equivalence is established.

## Output Tasks

- Generalize the compact `VarCorr()` output from the specialized `xar1`
  implementation.
- Do not print the full `member x time` correlation matrix by default.
- Provide full matrix output only for small blocks or explicit requests.
- Snapshot-test `print()`, `summary()`, `VarCorr()`, and
  `as.data.frame(VarCorr())`.

## Test Plan

Start from the specialized implementation as a reference branch.

1. Fixed-parameter objective equality:
   - specialized `homcsxar1()` vs generic separable `homcs x ar1`;
   - specialized `unxar1()` vs generic separable `us x ar1`.
2. Fit equality:
   - compare logLik, theta, fixed effects, fitted values, and `VarCorr()`.
3. Explicit covariance checks:
   - compare small-block full covariance matrices to
     `kronecker(R_time, Sigma_member)`.
4. Missing-data checks:
   - random dropout within groups;
   - whole member-time cells missing in some groups;
   - globally missing time point present in factor levels;
   - globally missing time point absent from levels should produce clear
     guidance.
5. Input validation checks:
   - incomplete grid;
   - mismatched margin/grid variables;
   - repeated coordinate variables;
   - intercept included accidentally;
   - unsupported margin combination;
   - both margins trying to carry scale.
6. Performance benchmarks:
   - specialized vs generic separable;
   - small dense Kronecker as a correctness-only comparison;
   - varying members, time points, groups, and missingness.

## Non-Goals For This Prototype

- Arbitrary `structure x structure` support.
- Ragged or group-specific grids.
- Irregular time.
- Role-specific AR parameters or VAR(1)-style nonseparable dynamics.
- Reduced-rank, proportional, or user-supplied covariance margins.
- Random slopes inside separable margins.

## PR Strategy

Keep the current specialized branch as a correctness/performance reference.
Develop this prototype from a clean upstream branch. If the generic backend can
match the specialized implementation in correctness and speed, the final PR
should contain the generic backend plus only the supported public structures.
