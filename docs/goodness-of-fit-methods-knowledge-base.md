# Goodness-of-fit and covariate adjustment methods inventory

Last reviewed: 2026-08-17

## Purpose and scope

This is a durable project knowledge base for evaluating and implementing
goodness-of-fit and covariate-adjusted effect sizes in WdStar. It records the
relevant behavior of the requested packages and related methods, based on
official documentation, package source, and primary papers.

The term “goodness of fit” is overloaded across these projects. This inventory
separates:

1. **Distance-variance fit:** sums of squares and R-squared computed from a
   distance/dissimilarity matrix.
2. **Ordination fit:** inertia explained by constrained axes, sites, or species.
3. **Predictive performance:** cross-validated error, AUC, Q-squared, or
   prediction R-squared.
4. **Differential-abundance or MANOVA inference:** test statistics and p-values
   that may adjust for covariates but do not necessarily report a model-fit
   measure.

The microbiome organization contains many packages. The review below covers the
components directly relevant to goodness-of-fit, beta-diversity, ordination,
PERMANOVA, and covariate adjustment: `microbiome`, `mia`, and the Orchestrating
Microbiome Analysis (OMA) workflow. It is not an inventory of every repository
under the organization.

## Terminology WdStar should use

Let:

- \(D\) be the original distance matrix;
- \(SS_T\) be its total distance-based sum of squares;
- \(Z\) be adjustment covariates;
- \(F\) be the focal group/model term;
- \(SS_E(Z)\) be residual sum of squares after fitting \(Z\); and
- \(SS_E(Z,F)\) be residual sum of squares after fitting \(Z+F\).

The following quantities answer different questions:

| Quantity | Definition | Interpretation |
|---|---:|---|
| Adjustment-only R-squared | \(R_Z^2=1-SS_E(Z)/SS_T\) | Fraction of original distance variation explained by adjustment covariates. |
| Full-model R-squared | \(R_{Z,F}^2=1-SS_E(Z,F)/SS_T\) | Fraction explained jointly by adjustment covariates and the focal term. |
| Semi-partial contribution of F after Z | \(R_{F\cdot Z,\,semi}^2=[SS_E(Z)-SS_E(Z,F)]/SS_T\) | Increment in full-model R-squared attributable uniquely to F after adjustment. It is on the original-total scale. |
| Partial R-squared of F given Z | \(R_{F\mid Z}^2=[SS_E(Z)-SS_E(Z,F)]/SS_E(Z)\) | Fraction of variation remaining after adjustment that F explains. |

The identity

\[
R_{F\cdot Z,\,semi}^2=R_{Z,F}^2-R_Z^2
\]

holds when all sums of squares are computed from the same original geometry and
the reduced and full models are genuinely nested.

Two uses of “adjusted R-squared” must not be conflated:

- **covariate-adjusted effect:** the contribution of F after controlling for Z;
- **degrees-of-freedom-adjusted R-squared:** a bias correction such as
  \(1-(1-R^2)(n-1)/(n-p-1)\).

`vegan::adonis2()` calls its term column “partial R2,” but its documentation
defines the values as term sum of squares divided by **total** sum of squares.
Under the terminology above, a conditional term value with that denominator is
a semi-partial contribution, not strict partial R-squared.

### Sequential versus incremental

An **incremental** effect is any specified contrast between two nested models.
For example, the increment of F after Z compares `Z + F` with `Z`. It describes
the estimand, not a general rule for assigning all terms.

A **sequential** analysis is a series of such increments in formula order. For
`Z1 + Z2 + F`, it compares intercept to `Z1`, `Z1` to `Z1 + Z2`, and then
`Z1 + Z2` to `Z1 + Z2 + F`. The results depend on order when predictors overlap.
Thus the final sequential term is exactly the desired increment of F after all
earlier Z terms, but earlier rows are not order-invariant marginal/Type III
effects. “Incremental” and “sequential” are therefore not opposites: sequential
analysis is one ordered way to generate multiple incremental contrasts.

## Package and method comparison

| Project | Primary purpose | Fit/effect-size implementation | Covariates included in that computation? | Relevance to WdStar |
|---|---|---|---|---|
| vegan | Community ecology, PERMANOVA, constrained ordination | `adonis2` distance R-squared; `dbrda`/`rda` inertia; `RsquareAdj`; `varpart` | **Yes.** Formula terms, `Condition()`, sequential and marginal tests, and partial variation partitions | Strongest open-source reference implementation |
| fast.adonis | Faster/weighted distance-based ANOVA | Projection-based weighted R-squared, sequential or marginal term SS, bootstrap SE | **Yes.** Covariates are ordinary RHS terms; sequential or leave-one-term-out marginal contributions | Useful source for direct matrix calculations and uncertainty, but not a full replacement for vegan |
| DISTLM / PERMANOVA+ | Distance-based multivariate regression and model selection | Direct Gower-matrix SS, R-squared, adjusted R-squared, AIC/AICc/BIC; marginal and conditional tests | **Yes.** Explicit reduced/full conditional tests; predictors may be fitted in sets | Closest conceptual precedent for nested covariate and focal models |
| MultANOVA | High-dimensional designed numeric responses | Multivariate Type III term R-squared from an effect orthogonalized to all other terms, divided by total SS | **Yes.** Each term is adjusted for all other formula terms | Strong precedent for order-invariant, total-denominator semi-partial contribution; not distance-based and design columns must be factors |
| mia / microbiome / OMA | Microbiome data containers and workflows | `mia::getPERMANOVA()` wraps `vegan::adonis2`; `getRDA()` wraps constrained ordination; legacy tutorials call vegan | **Yes, through vegan formulas.** `getPERMANOVA` defaults to marginal tests and also runs dispersion checks | Interface and reporting precedent, not an independent SS algorithm |
| ANCOM-BC / ANCOM-BC2 | Taxon-wise differential abundance with bias correction | Coefficients, SE, W statistics, p/q values, and global contrasts | **No comparable community-level R-squared.** Covariates are included in each taxon’s inferential model | Evidence that covariate adjustment belongs in the model, but not a GOF implementation to copy |
| coda4microbiome | Compositional predictive signatures | Apparent AUC/R-squared; CV AUC/MSE; survival C-index | **Yes, but predictively.** Covariates enter as an offset; reported fit combines offset and microbial signature | Not distance-based and does not isolate adjustment-only versus incremental focal R-squared |
| mixOmics | Latent-variable prediction and multi-omics integration | Explained component variance; CV error/BER/AUC; regression MSEP/RMSEP/R-squared/Q-squared | **Not as a general confounder-adjusted R-squared.** Multilevel decomposition handles paired/repeated designs, which is a different operation | Useful predictive-validation practices, not a template for distance effect-size partitioning |
| HDMedian | High-dimensional MANOVA via geometric medians and bootstrap | Max-type mean/median test statistics and bootstrap null distributions | **No.** Input is a list of group matrices; no covariate model or R-squared | Robust testing reference only |
| Nearing DA comparison repository | Benchmark of differential-abundance methods across datasets | FDR/FPR/consistency comparisons; an Aitchison PERMANOVA screening result | **No in the published pipeline.** Input metadata is described as sample ID plus one grouping, with at most two groups | Benchmark evidence, not a GOF algorithm |

## Implementation findings by project

### vegan

`vegan::adonis2()` accepts a distance matrix or community matrix on the left and
continuous variables, factors, transformations, and interactions on the right.
It offers:

- `by = "terms"`: sequential Type I contributions, first to last;
- `by = "margin"`: each term tested in a model containing all other terms;
- `Condition(Z)`: Z is partialled out before the constrained model is fit; and
- permutation tests under a reduced model.

Its R-squared column is `SumOfSqs / total SumOfSqs`. Therefore, for
`D ~ F + Condition(Z)`, the F value is a conditional, total-denominator
semi-partial contribution. A strict partial R-squared requires dividing the same
incremental SS by the residual SS of the Z-only model.

For constrained ordination, `RsquareAdj()` supplies both raw and
degrees-of-freedom-adjusted R-squared for RDA/dbRDA objects. `varpart()` accepts
dissimilarities and two to four explanatory tables, reports unique and shared
fractions, and primarily uses adjusted R-squared because raw variation fractions
are upward biased.

`vegan::goodness()` is not the whole-model statistic WdStar is seeking. It gives
cumulative inertia accounted for by selected axes for individual sites or
species. `inertcomp()` decomposes partial, constrained, and unconstrained inertia
at site/species level.

Sources:

- [`adonis2` reference](https://vegandevs.github.io/vegan/reference/adonis.html)
- [`RsquareAdj` reference](https://vegandevs.github.io/vegan/reference/RsquareAdj.html)
- [`varpart` reference](https://vegandevs.github.io/vegan/reference/varpart.html)
- [`goodness.cca` reference](https://vegandevs.github.io/vegan/reference/goodness.cca.html)
- [McArdle & Anderson (2001)](https://doi.org/10.1890/0012-9658(2001)082%5B0290:FMMTCD%5D2.0.CO;2)

### fast.adonis

The package expects \(A=-D^2/2\) rather than a `dist` object. Its source builds a
model matrix, handles aliased columns with QR, and computes weighted projection
R-squared. In simplified unweighted notation it uses the same centered-distance
idea as Gower-based distance regression.

For multiple terms, the source calculates:

- a full-model R-squared;
- marginal contributions as full R-squared minus a model omitting one term; and
- sequential contributions as differences between cumulative nested models.

The output converts those R-squared contributions back to sums of squares and
then reports `R2 = SumsOfSqs / total SS`. Thus its conditional term values also
use the total denominator. It can bootstrap the R-squared estimates and supports
survey/sample weights. Covariate adjustment is achieved by including covariates
as RHS terms and choosing their order for sequential analysis or using marginal
analysis.

Important limitations for adoption are the use of explicit `solve()` calls,
internal `vegan:::` functions, and the absence of a `Condition()` interface.

Sources:

- [fast.adonis repository](https://github.com/jennylsl/fast.adonis)
- [`fast.adonis.R` source](https://github.com/jennylsl/fast.adonis/blob/main/R/fast.adonis.R)
- [`suppli.code.R` source](https://github.com/jennylsl/fast.adonis/blob/main/R/suppli.code.R)

### DISTLM and PERMANOVA+

DISTLM is not an R package, but it is a primary methodological and software
precedent. Starting with the Gower-centered matrix \(G\) and a predictor hat
matrix \(H\), the manual defines:

\[
SS_T=\operatorname{tr}(G),\qquad
SS_{Reg}=\operatorname{tr}(HGH),\qquad
R^2=SS_{Reg}/SS_T.
\]

It distinguishes a marginal test of one predictor alone from a conditional test
of \(X_2\) after \(X_1\). The conditional pseudo-F explicitly compares a full
model \(X_1+X_2\) against a reduced model \(X_1\), and inference uses permutation
of residuals under the reduced model. The documentation explicitly describes
\(X_1\) as covariates in this setting.

DISTLM also supports predictor sets, user-specified order, forward/backward/
stepwise/best model search, and R-squared, adjusted R-squared, AIC, AICc, or BIC
selection criteria. These model-selection scores should not be confused with an
effect size for F after Z.

Sources:

- [PRIMER 8 DISTLM feature description](https://www.primer-e.com/features/primer-8-with-permanova/)
- [DISTLM partitioning and R-squared](https://learninghub.primer-e.com/books/permanova-for-primer-guide-to-software-and-statistical-methods/page/43-partitioning)
- [DISTLM conditional tests](https://learninghub.primer-e.com/books/permanova-for-primer-guide-to-software-and-statistical-methods/page/45-conditional-tests)
- [PERMANOVA designs with covariates](https://learninghub.primer-e.com/books/permanova-for-primer-guide-to-software-and-statistical-methods/page/135-designs-with-covariates-holdfast-invertebrates-revisited)

### MultANOVA

`MultANOVA()` fits a multivariate linear model to a numeric response matrix. For
each term, it constructs the term’s fitted effect and projects that effect away
from the design matrix containing every other term. If \(H_j\) is the resulting
Type III hypothesis cross-product matrix, its reported term effect is:

\[
R_j^2=\operatorname{tr}(H_j)/SS_T.
\]

This is order-invariant adjustment for the other formula terms and uses a total
denominator. It is therefore a clear implementation of a multivariate
semi-partial contribution, even though the package calls it “multivariate Type
III r-squared.”

The implementation is not directly reusable for WdStar because the response is
a rectangular numeric matrix rather than an arbitrary distance matrix, and the
source rejects non-factor design columns. Its projection logic is nevertheless
useful as a semantic reference.

Sources:

- [MultANOVA reference manual](https://stat.ethz.ch/CRAN/web/packages/MultANOVA/MultANOVA.pdf)
- [`MultANOVA.R` source](https://github.com/cran/MultANOVA/blob/master/R/MultANOVA.R)
- [Mahieu & Cariou (2025)](https://doi.org/10.1002/cem.70039)

### mia, microbiome, and OMA

`mia::getPERMANOVA()` delegates to `vegan::adonis2()`. It passes formula and
additional arguments through, defaults to `by = "margin"`, returns R-squared as
explained variance, and by default pairs the test with a `betadisper` homogeneity
assessment. Its data argument is explicitly described as including the
covariates defined by the formula.

The `microbiome` tutorials and OMA workflows use the same vegan/mia machinery.
Consequently, their support for covariate-adjusted distance fit is inherited
rather than an independent implementation.

Sources:

- [`mia::getPERMANOVA`](https://microbiome.github.io/mia/reference/getPERMANOVA.html)
- [microbiome PERMANOVA tutorial](https://microbiome.github.io/tutorials/PERMANOVA.html)
- [Orchestrating Microbiome Analysis](https://microbiome.github.io/OMA/)

### ANCOM-BC and ANCOM-BC2

ANCOM-BC is a taxon-wise differential-abundance framework, not a global
distance-variance model. `ancombc()` uses a formula describing abundance as a
function of metadata. `ancombc2()` uses `fix_formula` and optionally random
effects. The output includes log-fold changes, standard errors, test statistics,
p-values, q-values, detection flags, and global/multiple-comparison results.

The methods therefore **do adjust inference for covariates**, but the reviewed
public API does not return a community-level distance R-squared, an adjustment-
only R-squared, or an incremental group R-squared. It is not a goodness-of-fit
implementation for WdStar to emulate.

Sources:

- [ANCOMBC Bioconductor manual](https://bioconductor.org/packages/release/bioc/manuals/ANCOMBC/man/ANCOMBC.pdf)
- [ANCOMBC source repository](https://github.com/FrederickHuangLin/ANCOMBC)

### coda4microbiome

`coda_glmnet()` fits elastic-net models to pairwise log ratios. For binary
outcomes it reports apparent and cross-validated AUC. For continuous outcomes it
reports apparent `cor(predictions, y)^2` and cross-validated MSE. `coda_coxnet()`
reports apparent and cross-validated C-index.

When `covar` is supplied, the source first fits an ordinary GLM/LM to the
covariates, then supplies its fitted linear predictor as an offset to
`cv.glmnet()`. Final predictions include that offset. Therefore the apparent
R-squared or AUC is for the **combined covariate plus microbial-signature
prediction**, not the microbial signature’s unique contribution after
adjustment. The package does not report a Z-only fit and does not compute either
of WdStar’s nested-model incremental ratios.

Sources:

- [`coda_glmnet` reference](https://search.r-project.org/CRAN/refmans/coda4microbiome/html/coda_glmnet.html)
- [`coda_coxnet` reference](https://search.r-project.org/CRAN/refmans/coda4microbiome/html/coda_coxnet.html)
- [`coda_glmnet` source](https://github.com/cran/coda4microbiome/blob/master/R/coda4microbiome_functions.R)
- [Calle et al. (2023)](https://doi.org/10.1186/s12859-023-05171-4)

### mixOmics

mixOmics reports multiple quantities that should not be treated as distance-
based goodness of fit:

- component-wise variance explained for PCA and latent-variable models;
- classification error rate, balanced error rate, and AUC;
- MSEP/RMSEP, R-squared, and Q-squared for PLS regression; and
- cross-validation from `perf.assess()`/`perf()` and tuning functions.

The general PLS/PLS-DA interfaces do not provide a formula/`Condition()`-style
confounder adjustment whose R-squared is partitioned into Z-only and focal
incremental components. The multilevel option removes between-subject variation
for repeated or paired measurements; that is repeated-measures decomposition,
not general covariate-adjusted goodness of fit.

Sources:

- [mixOmics performance assessment](https://mixomics.org/performance-assessment-and-parameter-tuning/)
- [mixOmics PLS methods](https://mixomics.org/methods/spls/)
- [mixOmics multilevel methods](https://mixomics.org/methods/multilevel/)
- [mixOmics source repository](https://github.com/mixOmics-org/mixOmics)

### HDMedian

The exposed `MANOVA_Median()` implementation receives a list of group-specific
matrices. It computes geometric medians, maximum pairwise standardized or
unstandardized median differences, and multiplier-bootstrap replicates. Related
functions do the same for means.

There is no model matrix, no covariate argument, no reduced/full model
comparison, and no R-squared. It contributes a robust high-dimensional test,
not a covariate-adjusted GOF method.

Sources:

- [HDMedian repository and implementation](https://github.com/liuhuapeng666/HDMedian)
- [Liu et al. (2024)](https://doi.org/10.1093/biomtc/ujae088)

### Nearing comparison repository

This repository is an analysis pipeline for comparing the outputs of multiple
differential-abundance methods across 38 datasets. Its evaluation targets method
agreement, false-positive behavior, and consistency. It stores Aitchison-
distance PERMANOVA results as a dataset-level screening analysis, but it does not
define a reusable GOF function.

The documented pipeline expects two-column metadata containing sample IDs and a
single grouping with at most two groups. It therefore does not provide a
covariate-adjusted goodness-of-fit implementation.

Sources:

- [Comparison repository](https://github.com/nearinj/Comparison_of_DA_microbiome_methods)
- [Nearing et al. (2022)](https://doi.org/10.1038/s41467-022-28034-z)

## Conclusions for WdStar

### Which reviewed methods adjust goodness of fit for covariates?

**Directly comparable distance or multivariate variance partitioning:**

- `vegan` (`adonis2`, `dbrda`, `varpart`);
- `fast.adonis` through formula terms and nested/marginal fits;
- DISTLM through reduced/full conditional tests;
- `MultANOVA` through Type III adjustment for all other terms; and
- `mia`/`microbiome`/OMA through their use of vegan.

**Covariates are modeled, but the reported quantity is not comparable to WdStar
distance R-squared:**

- ANCOM-BC/ANCOM-BC2: covariate-adjusted taxon-wise inference, no global GOF;
- coda4microbiome: combined predictive fit with covariate offset; and
- mixOmics: predictive/latent-variable fit and repeated-measures decomposition,
  not general confounder-partitioned distance R-squared.

**No relevant covariate-adjusted GOF implementation:** HDMedian and the Nearing
benchmark repository.

### Recommended computational core

For a WdStar implementation, use the direct distance geometry employed by
McArdle-Anderson, vegan, and DISTLM:

1. Validate and align the original distance matrix and model data once.
2. Form \(A=-D^2/2\), \(J=I-11'/n\), and \(G=JAJ\) once.
3. Build full-rank design matrices for the intercept-only, Z-only, and Z+F
   models using pivoted QR. Avoid explicit matrix inversion.
4. Compute projection or residual operators from the QR bases and obtain
   \(SS_T\), \(SS_E(Z)\), and \(SS_E(Z,F)\) from the same Gower matrix.
5. Return, with unambiguous names, adjustment-only R-squared, full-model
   R-squared, the semi-partial contribution of F, and strict partial R-squared.
6. If a degrees-of-freedom-adjusted R-squared is offered, label it separately
   and report the rank/df used.
7. Treat permutation restrictions (`strata`) separately from adjustment
   covariates. For a covariate-adjusted significance test of F, use a
   reduced-model residual permutation scheme such as Freedman-Lane; do not
   permute raw labels across a covariate structure.
8. Report a dispersion diagnostic alongside centroid/location inference where
   appropriate, following vegan/mia practice.

### Important implication for the current `a.dist()` route

The current `a.dist()` constructs a residual Gower-type matrix, eigen-decomposes
it, and sets negative eigenvalues to zero before reconstructing Euclidean
distances. Repeating that operation independently for Z and Z+F can change the
geometry differently in the two fits. Consequently,

```r
dist.goodness.of.fit(D, a.dist(D, ~ Z + F)) -
  dist.goodness.of.fit(D, a.dist(D, ~ Z))
```

is a useful provisional semi-partial calculation, but it is not guaranteed to
equal an exact nested-model SS partition for non-Euclidean dissimilarities.
The production implementation should compute both fits against one shared Gower
matrix. If Euclidification is chosen, apply and document it once at the original
distance level rather than truncating each residual model independently.

### Future improvement: pre-adjusted distances

`WdS.test(dm = a.dm, f = ...)` currently preserves any
`distance.diagnostics` attribute attached to a pre-adjusted distance matrix as
role `"input"`. This records how the supplied distance matrix was constructed,
but it does not by itself provide enough information to recover goodness-of-fit
values relative to the original, unadjusted distance matrix.

A low-memory future extension could have `a.dist()` store scalar variation
metadata alongside the existing eigendecomposition diagnostics, for example the
original distance variation, the residual distance variation, the adjustment
formula, and the tolerance used during reconstruction. With those scalars,
`WdS.test(dm = a.dm, f = ..., goodness = "adjustment")` could recover the
adjustment-only pseudo-\(R^2\),

```text
1 - residual_variation / original_variation
```

without storing the original distance matrix inside `a.dm`.

This would still not be enough to compute factor-only, full-model,
semi-partial, or partial goodness-of-fit for `f`. Those components require the
original distance matrix so the factor-only and combined adjustment-plus-factor
residual distances can be constructed against the same original geometry.
Supporting those values for externally adjusted inputs would require a more
explicit API, such as accepting both the original distance matrix and the
pre-adjusted distance matrix. The package should avoid storing the full original
distance matrix as a default attribute because that would make simulation and
large-data workflows unnecessarily memory-heavy.

### Naming recommendation

Avoid a bare `partial.R2` field. Prefer explicit output names such as:

```text
adjustment.r.squared
full.model.r.squared
group.semi.partial.r.squared
group.partial.r.squared
df.adjusted.r.squared       # only if implemented
```

This prevents the terminology mismatch found in existing software from becoming
part of WdStar’s public API.
