

# $W_d^*$: distance-based multivariate analysis of variance for multivariate data

## Introduction 

`WdStar` is an R package for $W_d^*$, a method for multivariate analysis of variance based on Welch's MANOVA designed to address challenges of existing methods including PERMANOVA.
- Are you using PERMANOVA and concerned about how heteroscedasticity and unbalanced sample sizes may affect your analyses? 
- Have you ever questioned the reliability of community-wide microbiome analysis results?
- Are you looking for a global test for your multivariate dataset?

>**$W_d^*$ addresses these concerns and is**
>- robust to heteroscedasticity; 
>- handles multi-level factors and stratification;
>- allows for multiple post hoc testing scenarios;
>- allows for adjustment of covariates; and
>- can be used with any distance or dissimilarity matrix.


## Peer-Reviewed Publications on the $W_d^*$-test Family  
**Preprint: Covariate-adjusted multivariate analysis for omics data**  
- Hamidi B, Fanning L, Wallace K, & Alekseyenko AV. *Submitted to Briefings in Bioinformatics on June 23, 2026. Not currently accepted or under review.* DOI forthcoming.
- [Code repository](https://github.com/alekseyenko/WdStar/tree/master/publications/Hamidi%20et%20al.%20Bioinformatics%20(2026))


**$W_d^*$-test: robust distance-based multivariate analysis of variance**  
- [Hamidi B, Wallace K, Vasu C, & Alekseyenko AV. *Microbiome.* 2019.](https://doi.org/10.1186/s40168-019-0659-9)
- [Code repository](https://github.com/alekseyenko/WdStar/tree/master/publications/Hamidi%20et%20al.%20Microbiome%20(2019))


**$T_w^2$: multivariate Welch t-test on distances**  
- [Alekseyenko AV. Multivariate Welch t-test on distances. *Bioinformatics*. 2016.](https://doi.org/10.1093/bioinformatics/btw524) 
- [Code repository](https://github.com/alekseyenko/Tw2)

## Methods Knowledge Base

See the [goodness-of-fit and covariate-adjustment methods inventory](https://github.com/alekseyenko/WdStar/blob/main/docs/goodness-of-fit-methods-knowledge-base.md) for implementation notes and comparisons with vegan, fast.adonis, DISTLM, MultANOVA, and microbiome-analysis packages.

## Installation  
Source installation of the `WdStar` R package is available directly from GitHub using `remotes` for R 3.6 or later:
```R
install.packages("remotes")
remotes::install_github("alekseyenko/WdStar", force = TRUE)
library(WdStar)
packageVersion("WdStar")
```

Until a CRAN release is available, install `WdStar` from GitHub as shown above.


## Quick Start  

For detailed and complex examples please refer to [our publication repositories](https://github.com/alekseyenko?tab=readme-ov-file#peer-reviewed-publications-on-the-w_d-test-family), which contain Markdown files with application datasets and code.


The following is a simple example using the `mtcars` dataset to assess the effect of `gear` on `mpg`, `cyl`, and `disp` (first three variables of the dataset):   

```R
# Load dataset
data(mtcars)

# The outcome could be a single variable or multiple variables (such as multidimensional omics data).  

### This is an example of outcome with a single variable (`mpg`):
dm <- dist(mtcars$mpg, method="euclidean")

### This is an example of outcome with multiple variables (`mpg`, `cyl`, and `disp`):
dm <- dist(mtcars[1:3], method="euclidean") 

# Grouping/independent variable. You could use multiple variables here too.
f <- factor(mtcars$gear)

# Basic multivariate test example ###########
#############################################
WdS.test(dm=dm, f=f)

## By default, the unadjusted test's goodness.of.fit reports the
## distance-based pseudo-R-squared for the tested factor.
unadjusted_res <- WdS.test(dm=dm, f=f)
unadjusted_res$goodness.of.fit

## Use goodness="none" to skip goodness-of-fit calculations.
WdS.test(dm=dm, f=f, goodness="none")

# Stratified example ########################
#############################################
strata <- factor(mtcars$vs)
WdS.test(dm=dm, f=f, strata=strata)

# Covariate adjustment/elimination examples #
#############################################
## Right-hand side adjustment formula to specify adjustment covariates. 
formula <- ~ wt + as.factor(am) 

## Adjustment example 1: pass unadjusted `dm` and formula to WdS.test()
WdS.test(dm=dm, f=f, formula=formula, formula_data=mtcars) ## Perform adjusted test

## By default, goodness.of.fit reports the distance-based semi-partial
## pseudo-R-squared for the tested factor after adjustment.
## Additional components computed along the way are stored but not printed.
res <- WdS.test(dm=dm, f=f, formula=formula, formula_data=mtcars)
res$goodness.components

## Eigenvalue and tolerance diagnostics for residual distance matrices are
## stored separately from goodness-of-fit values.
res$distance.diagnostics

## Request all available goodness-of-fit components, including factor-only,
## adjustment-only, full-model, semi-partial, and partial pseudo-R-squared.
WdS.test(dm=dm, f=f, formula=formula, formula_data=mtcars, goodness="all")

## Interpretation of adjusted components:
## - adjustment: variation explained by the adjustment variables alone.
## - full: variation explained jointly by adjustment variables and the tested factor.
## - semi.partial: additional total variation explained by the tested factor after adjustment.
## - partial: adjustment-residual variation explained by the tested factor.
## Negative values can occur if the residual distances contain more variation
## than the original distance matrix.

## Adjustment example 2: Create the adjusted distance matrix `a.dm` outside the function
a.dm <- a.dist(dm=dm, formula=formula, formula_data=mtcars) 
WdS.test(dm=a.dm, f=f) ## Perform adjusted test with `a.dm`; input diagnostics are preserved.
attr(a.dm, "distance.diagnostics")

## Store raw eigenvalues only when deeper diagnostics are needed.
a.dm.with.eigenvalues <- a.dist(
  dm=dm,
  formula=formula,
  formula_data=mtcars,
  keep.eigenvalues=TRUE
)
length(attr(a.dm.with.eigenvalues, "distance.diagnostics")$eigenvalues[[1]])

## Distance-based pseudo-R-squared can also be computed directly
dist.goodness.of.fit(dm=dm, dm_residual=a.dm)
```

Further examples are provided in the package documentation and may be accessed by running the following commands:
```R
?WdS.test
?a.dist
?dist.goodness.of.fit
```

## Feature Requests and Bugs
We welcome feature requests and bug reports and kindly ask you to submit them via [our GitHub issue tracker.](https://github.com/alekseyenko/WdStar/issues)
