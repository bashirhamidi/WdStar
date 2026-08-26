# WdStar 2.4.0 2026-08-26

* Added `dist.goodness.of.fit()` to compute a distance-based pseudo-coefficient
  of determination, \(R^2_{pseudo} = 1 - \sigma^2_{residual} /
  \sigma^2_{total}\), for users who compute original and residual distance
  matrices outside of `WdS.test()`. The helper now uses the more general
  `dm_residual` argument, while keeping `adjusted_dm` as a compatibility alias,
  and validates that both inputs are `dist` objects with the same number of
  observations.

* `WdS.test()` now reports a separate `goodness.of.fit` component. This keeps
  `estimate` reserved for omega-squared (\(\omega^2\)) while reporting
  distance-based pseudo-\(R^2\) separately. The new `goodness` argument defaults
  to `"auto"`, which reports the tested factor's pseudo-\(R^2\) for unadjusted
  tests and the semi-partial pseudo-\(R^2\) for adjusted tests.

* Added `goodness.components` to `WdS.test()` results. This table stores the
  goodness-of-fit components computed along the way, including adjustment-only,
  full-model, semi-partial, and partial pseudo-\(R^2\) for adjusted tests.
  Users can request specific components or all available components with
  `goodness = "all"`.

* Added compact distance diagnostics to `a.dist()` and `WdS.test()`.
  `a.dist()` now attaches a `distance.diagnostics` attribute that records the
  eigendecomposition tolerance, raw positive/negative eigenvalue counts,
  eigenvalues zeroed by `tol`, and negative eigenvalues removed during distance
  reconstruction. `WdS.test()` collects these diagnostics in a separate
  `distance.diagnostics` output component for constructed residual distance
  matrices. Raw eigenvalue vectors are stored only when
  `keep.eigenvalues = TRUE`.

* Improved eigenvalue handling in `a.dist()`: eigenvalues smaller than `tol` and
  negative eigenvalues are now reported to the user and set to zero before the
  adjusted distance matrix is reconstructed.

* Improved `formula_data` handling for `a.dist()` and adjusted `WdS.test()`
  calls. Formula variables can be resolved from the caller environment by
  default, and `formula_data` now accepts environments, data frames, lists, and
  objects coercible to data frames, including `phyloseq::sample_data()` objects.

* Added `phyloseq` as a suggested package to support tests and examples involving
  `sample_data()` inputs.

* Updated documentation for `a.dist()`, `WdS.test()`, and
  `dist.goodness.of.fit()`, including accepted `formula_data` input types,
  parent-frame formula lookup, goodness-of-fit component definitions, distance
  diagnostics, and direct pseudo-\(R^2\) calculation from original and residual
  distance matrices.

* Added regression tests for parent-frame formula lookup,
  `phyloseq::sample_data()` compatibility, adjusted-test goodness-of-fit, direct
  `dist.goodness.of.fit()` use, `a.dist()` tolerance messages, and negative
  eigenvalue handling.

* Added GitHub issue templates for bug reports and feature requests.

# WdStar 2.3.0 2025-09-21

* Added package versioning.

* Documented installation from GitHub using `remotes::install_github()`.

* Added covariate-adjusted testing support through `WdS.test()`.

* Added omega-squared as an effect size estimate.

* Added between-group degrees of freedom to `WdS.test()` output.

* Updated `a.dist()` parameters to avoid requesting duplicate objects.

* Improved error handling and data type handling.

* Improved `a.dist()` and `WdS.test()` outputs.

* Revamped help files and updated examples for `a.dist()` and `WdS.test()`.
