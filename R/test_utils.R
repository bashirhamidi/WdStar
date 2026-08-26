#' Conducts a generic distance-based permutation test for k-group differences
#'
#' This function performs a distance-based permutation test using an arbitrary
#' test statistic. It is built on the principle of Welch's ANOVA and extends it
#' to handle multivariate distance data. It is particularly useful for analyzing
#' microbiome data.
#'
#' @param test.statistic A function to calculate the test statistic.
#' @param dm A distance metric (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#' @param strata A factor variable representing the strata. Default is NULL.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{statistic}: The observed test statistic
#'   \item \code{nrep}: The number of permutations performed
#' }
#'
#' @references Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based multivariate analysis of variance." Microbiome 7.1 (2019):1-9.
#' @seealso \code{\link{Tw2.test}}, \code{\link{WdS.test}}
#' @export
#' @examples
#' # TODO: Add examples
generic.distance.permutation.test =
  function(test.statistic, dm, f, nrep=999, strata = NULL){
    N = length(f)
    generate.permutation=function(){
      f[sample(N)]
    }

    if(!is.null(strata)){
      # map elements of each strata back to their positions in the factor variable
      strata.map = order(unlist(tapply(seq_along(f), strata, identity)))
      generate.permutation=function(){
        p = unlist(tapply(f,strata,sample)) # permute within strata
        p[strata.map]
      }
    }

    stats = c(test.statistic(dm, f),
              replicate(nrep,
                        test.statistic(dm, generate.permutation())))

    p.value = sum(stats>=stats[1])/(nrep+1)
    statistic = stats[1]
    list(p.value = p.value, statistic = statistic, nrep=nrep)
  }

#' Conducts a Tw2 distance-based permutation test for k-group differences
#'
#' This function is a specialized version of the \code{\link{generic.distance.permutation.test}},
#' specifically designed to use the Tw2 test statistic for k-group comparison.
#'
#' @param dm A distance metric (any arbitrary distance or dissimilarity metric)
#' @param f A factor variable representing the groups.
#' @param nrep The number of permutations to conduct. Default is 999.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{statistic}: The observed test statistic
#'   \item \code{nrep}: The number of permutations performed
#' }
#'
#' @seealso \code{\link{WdS.test}}, \code{\link{generic.distance.permutation.test}}
#' @export
#' @examples
#' # TODO: Add examples
Tw2.test <- function(dm, f, nrep = 999) {
  generic.distance.permutation.test(Tw2, dm = dm, f = f, nrep = nrep)
}

#' Conducts distance-based multivariate Welch ANOVA
#'
#' This function performs the \eqn{\mathnormal{W}_d^*} test statistic for
#' k-group comparison using a given distance matrix.
#'
#' This method uses permutation testing to establish the significance by computing
#' \eqn{\mathnormal{W}_d^*} statistic on \eqn{\mathnormal{m}} permutations of the original data
#' and estimate the significance as the fraction of times the permuted statistic
#' is greater than or equal to \eqn{\mathnormal{W}_d}.
#'
#' Our approach can accommodate multi-level factors, stratification (via
#' restricted permutations), multiple post-hoc testing scenarios, and covariate
#' adjustment/ elimination (via projection of residuals).
#'
#' If optional covariate adjustment parameters (\code{formula}, and
#' \code{formula_data}) are provided as input, they are passed to the \code{WdStar::a.dist()}
#' function to project residual matrices and remove the effect of specified
#' covariate(s) from distance matrix \code{dm} before performing the test.
#' Users may choose to make adjustments within the function or use \code{WdStar::a.dist()}
#' outside the function and then pass on the adjusted distance matrix to \code{WdStar::WdS.test()}.
#'
#' Goodness-of-fit is reported as a distance-based pseudo-\eqn{R^2}. With
#' \code{goodness = "auto"}, unadjusted tests report \code{factor.only}, the
#' variation explained by the tested factor. Adjusted tests report
#' \code{semi.partial}, the additional fraction of total distance variation
#' explained by the tested factor after the adjustment variables. The
#' \code{goodness.components} table stores components computed directly or
#' derived from the residual distance matrices used for the requested
#' calculation. For adjusted \code{"auto"} tests, this includes
#' \code{adjustment}, \code{full}, \code{semi.partial}, and \code{partial}, but
#' not \code{factor.only}, which requires an additional residual distance matrix.
#' Values near zero indicate little distance-variation reduction, larger
#' positive values indicate greater reduction, and negative values may occur
#' when the residual distance matrix has more variation than the original
#' distance matrix.
#'
#' Distance diagnostics are stored separately from goodness-of-fit. When
#' \code{distance.diagnostics = TRUE}, \code{WdS.test()} collects the diagnostic
#' summaries attached by \code{\link{a.dist}()} for each residual distance
#' matrix it constructs, using roles such as \code{"adjustment"},
#' \code{"factor.only"}, and \code{"full"}. These diagnostics record the
#' tolerance used, how many eigenvalues were treated as zero by \code{tol}, and
#' how many negative eigenvalues were removed during reconstruction. Use
#' \code{keep.eigenvalues = TRUE} to store the raw eigenvalue vectors for deeper
#' diagnosis; the default keeps only the compact summary.
#'
#' The diagnostic rows correspond to residual distance matrices actually
#' constructed by the call. For example, an adjusted test with
#' \code{goodness = "auto"} stores rows for \code{"adjustment"} and
#' \code{"full"}; \code{"factor.only"} is added only when requested, such as
#' with \code{goodness = "all"}. An unadjusted test with the default
#' \code{goodness = "auto"} stores a \code{"factor.only"} row. If
#' \code{goodness = "none"} and no adjustment formula is supplied, no residual
#' distance matrix is constructed and \code{distance.diagnostics} is
#' \code{NULL}, unless the input \code{dm} already carries diagnostics from a
#' previous \code{a.dist()} call. In that case, the input diagnostics are
#' preserved with role \code{"input"}.
#'
#' For adjustment variables \eqn{Z} and tested factor \eqn{F}, the components
#' are interpreted as follows:
#' \itemize{
#'   \item \code{factor.only}: \eqn{R^2_F}, variation explained by \eqn{F}
#'   alone, ignoring adjustment variables.
#'   \item \code{adjustment}: \eqn{R^2_Z}, variation explained by \eqn{Z}
#'   alone, ignoring \eqn{F}.
#'   \item \code{full}: \eqn{R^2_{Z+F}}, variation explained by \eqn{Z} and
#'   \eqn{F} together from one model formula.
#'   \item \code{semi.partial}: \eqn{R^2_{Z+F} - R^2_Z}, the additional
#'   fraction of total distance variation explained by \eqn{F} after \eqn{Z}.
#'   \item \code{partial}: \eqn{(R^2_{Z+F} - R^2_Z)/(1 - R^2_Z)}, the fraction
#'   of adjustment-residual variation explained by \eqn{F}.
#' }
#' The full-model residual matrix is computed from the original \code{dm} with
#' a combined formula. It is not computed by applying \code{a.dist()} a second
#' time to an already residualized distance matrix.
#' When no adjustment formula is supplied, \code{adjustment} is unavailable; if
#' requested, \code{full}, \code{semi.partial}, and \code{partial} reduce to the
#' same value as \code{factor.only}. Use \code{goodness = "none"} to skip these
#' additional goodness-of-fit calculations.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param nrep (optional) The number of permutations to conduct. Default is 999.
#' @param strata (optional) A factor variable representing strata. If specified,
#'   the test will perform p-value computations using stratified permutation.
#'   Default is NULL.
#' @param formula (optional) A right-hand side ONLY formula to be used with
#'   \code{formula_data} for confounder adjustment/elimination from \code{dm}.
#'   Results in
#'   adjusted WdS statistic, omega-squared (\eqn{\omega^2}{\omega^2}) effect
#'   size estimate, distance-based pseudo-\eqn{R^2},
#'   and p-value. Note that you may use any data type including factor,
#'   character, integer, and numeric. Default is NULL
#' @param formula_data (optional) An environment, data frame, list, or object
#'   coercible to a data frame, such as \code{phyloseq::sample_data()}, to be
#'   used in conjunction with \code{formula} for confounder adjustment. Default is
#'   parent.frame()
#' @param goodness Which distance-based pseudo-\eqn{R^2} goodness-of-fit
#'   component to report. Options are \code{"auto"}, \code{"factor.only"},
#'   \code{"adjustment"}, \code{"full"}, \code{"semi.partial"},
#'   \code{"partial"}, \code{"all"}, and \code{"none"}. The default
#'   \code{"auto"} reports \code{"factor.only"} for unadjusted tests and
#'   \code{"semi.partial"} for adjusted tests. Multiple explicit components may
#'   be requested with a character vector, except for \code{"auto"},
#'   \code{"all"}, and \code{"none"}. Use \code{"none"} to skip
#'   goodness-of-fit calculations.
#' @param distance.diagnostics Logical indicating whether to store compact
#'   eigendecomposition diagnostics for residual distance matrices constructed
#'   inside \code{WdS.test()}. Default is \code{TRUE}.
#' @param keep.eigenvalues Logical indicating whether to include raw eigenvalue
#'   vectors in \code{distance.diagnostics}. Default is \code{FALSE} to avoid
#'   extra memory use in large analyses and simulations.
#' @param ... Additional arguments passed to \code{\link{a.dist}()}, such as
#'   \code{tol} (tolerance for eigenvalues; default \code{1e-8}). These
#'   arguments are used when residual distance matrices are computed for
#'   adjustment or goodness-of-fit components.
#' @return A list containing:
#' \itemize{
#'   \item \code{method}: name of the method used.
#'   \item \code{data.name}: a string describing the input data.
#'   \item \code{statistic}: observed WdS test statistic.
#'   \item \code{estimate}: effect size estimator of variance, omega-squared (\eqn{\omega^2}{\omega^2}).
#'   \item \code{goodness.of.fit}: requested/default distance-based
#'   pseudo-\eqn{R^2} component.
#'   \item \code{goodness.components}: data frame of the goodness-of-fit
#'   components computed while evaluating \code{goodness.of.fit}. This may
#'   include intermediate values that are not printed by default.
#'   \item \code{distance.diagnostics}: data frame of eigendecomposition
#'   diagnostics for residual distance matrices constructed by \code{a.dist()},
#'   or \code{NULL} when no diagnostics were collected.
#'   \item \code{p.value}: The p-value of the test.
#'   \item \code{parameter}: A list of 2-4 containing:
#'          \itemize{
#'          \item \code{strata}: (optional) strata variable used to perform restricted permutations.
#'          \item \code{formula}: (optional) formula for residual projection to perform matrix adjustment.
#'          \item \code{dfb}: between degrees of freedom.
#'          \item \code{nrep}: number of permutations performed
#'          }
#' }
#'
#' @seealso \code{\link{generic.distance.permutation.test}}, \code{\link{Tw2.test}},
#'          \code{\link{a.dist}}, \code{\link{dist.goodness.of.fit}}
#' @importFrom stats terms
#' @export
#' @examples
#' # The following is a simple example using the mtcars dataset to assess the
#' # effect of gears on mpg, cyl, and disp (first three variables of the dataset):
#' data(mtcars)
#'
#' # The outcome could be a single variable or multiple variables (such as
#' # multidimensional omics data).
#'
#' ## This is an example of outcome with a single variable:
#' dm <- dist(mtcars$mpg, method="euclidean")
#'
#' ## This is an example of outcome with multiple variables:
#' dm <- dist(mtcars[1:3], method="euclidean")
#'
#' # Grouping/independent variable. You could use multiple variables here too.
#' f <- factor(mtcars$gear)
#'
#' # Basic multivariate test example ###########
#' #############################################
#' WdS.test(dm=dm, f=f)
#'
#' # Stratified example ########################
#' #############################################
#' strata <- factor(mtcars$vs)
#' WdS.test(dm=dm, f=f, strata=strata)
#'
#' # Covariate adjustment/elimination examples #
#' #############################################
#' ## Right-hand side adjustment formula to specify adjustment covariates.
#' formula <- ~ wt + as.factor(am)
#'
#' ## Adjustment example 1: pass unadjusted `dm` and formula to WdS.test()
#' WdS.test(dm=dm, f=f, formula=formula, formula_data=mtcars) ## Perform adjusted test
#'
#' ## Inspect eigenvalue/tolerance diagnostics for constructed distance matrices
#' res <- WdS.test(dm=dm, f=f, formula=formula, formula_data=mtcars)
#' res$distance.diagnostics
#'
#' ## Request all available goodness-of-fit components
#' WdS.test(dm=dm, f=f, formula=formula, formula_data=mtcars, goodness="all")
#'
#' ## Adjustment example 2: Create the adjusted distance matrix `a.dm` outside
#' ## the function
#' a.dm <- a.dist(dm=dm, formula=formula, formula_data=mtcars)
#' WdS.test(dm=a.dm, f=f) ## Perform adjusted test with `a.dm`
#'
# goodness and diagnostic controls are placed after ... so older calls that pass
# a.dist() arguments through ... keep their positional behavior; users should
# pass these newer controls by name.
WdS.test <- function(dm, f, nrep=999, strata=NULL, formula=NULL, formula_data=parent.frame(), ...,
                     goodness="auto", distance.diagnostics = TRUE,
                     keep.eigenvalues = FALSE){
  .validate_distance_diagnostic_options(distance.diagnostics, keep.eigenvalues)
  has_formula <- !is.null(formula)
  has_formula_data <- !missing(formula_data) && !is.null(formula_data)

  # Goodness-of-fit always compares a residual distance matrix back to the
  # original dm, even when the test statistic itself uses an adjusted dm.
  original_dm <- dm
  data <- NULL

  if (!has_formula && has_formula_data) {
      stop("Both 'formula' and 'formula_data' must be provided together for a.dist() adjustment processing.")
   }
  if (has_formula) {
    data <- .as_formula_data(formula_data)
    dm <- a.dist(
      dm,
      formula,
      data,
      ...,
      distance.diagnostics = distance.diagnostics,
      keep.eigenvalues = keep.eigenvalues
    )
  }
  test.results <- generic.distance.permutation.test(WdS, dm=dm, f=f, nrep=nrep, strata=strata)

    # Name of the hypothesis test
  if(!is.null(formula)) {
    method <- "Adjusted Distance-based Multivariate Welch ANOVA"
    adjustment_formula <- deparse(substitute(formula))
    } else {
    method <- "Distance-based Multivariate Welch ANOVA"
    adjustment_formula <- NULL
    }

  # Strata information
  if(!is.null(strata)){
    strata_selected <- "\n\tP-value is based on stratified permutations. WdS statistic and \u03C9\u00B2 are not computed by strata."}
    else {strata_selected <- NULL}

  # Description of the data
  data.name <- paste0(attr(dm, "method"), " distance matrix with ", attr(dm, "Size"), " observations and grouping factor with ", nlevels(f), " levels", strata_selected)

  # Value of the parameter under the null hypothesis
  # null.value <- 0
  # attr(null.value, "names") <- "expected WdS under H0"

  # Direction of the alternative hypothesis relative to the null value
  # alternative <- "two.sided"

  # Statistic value
  statistic <- test.results$statistic
  attr(statistic, "names") <- "WdS"

  # Effect size estimator of variance, omega-squared
  estimate <- (((length(levels(f))-1)*(unname(statistic)-1))/((length(levels(f))-1)*(unname(statistic) - 1) + attr(dm, "Size")))
  attr(estimate, "names") <- "effect size estimator of variance, omega-squared (\u03C9\u00B2)"

  # Build the requested pseudo-R2 output after the WdS statistic is computed.
  # For adjusted tests, dm is now the adjustment-residual distance matrix.
  goodness <- .compute_wdstest_goodness(
    dm = original_dm,
    f = f,
    formula = formula,
    formula_data = data,
    dm_adjustment_residual = if (has_formula) dm else NULL,
    goodness = goodness,
    distance.diagnostics = distance.diagnostics,
    keep.eigenvalues = keep.eigenvalues,
    ...
  )

  # P-value
  p.value <- test.results$p.value

  # Creating object of class 'htest'
  TEST <- list(method = method, data.name = data.name, # null.value = null.value, # alternative = alternative,
               statistic = statistic, p.value = p.value,
               estimate = estimate,
               goodness.of.fit = goodness$goodness.of.fit,
               goodness.components = goodness$goodness.components,
               distance.diagnostics = goodness$distance.diagnostics,
               parameter = c(if(!is.null(strata)){c("strata"= deparse(substitute(strata)))},
                             if(!is.null(formula)){c("formula"= paste0(deparse(substitute(formula)),
                                                                       " with ",
                                                                       length(attr(terms(formula), "term.labels")),
                                                                       " terms"))},
                             c("dfb" = (length(levels(f))-1)),
                             c("number of permutations" = nrep)
                             )
               )

  class(TEST) <- c("wdstest", "htest")
  return(TEST)
}

#' @export
print.wdstest <- function(x, ...) {
  estimate_name <- .format_estimate_name(names(x$estimate))

  x$estimate <- stats::setNames(unname(x$estimate), estimate_name)
  NextMethod()

  # htest printing does not know about WdStar's goodness fields, so print them
  # after the standard htest block while keeping estimate reserved for omega2.
  if (!is.null(x$goodness.of.fit)) {
    goodness_names <- .format_r2_name(names(x$goodness.of.fit))
    if (length(x$goodness.of.fit) == 1) {
      cat(goodness_names, ": ", format(unname(x$goodness.of.fit), digits = 4), "\n", sep = "")
    } else {
      cat("goodness-of-fit:\n")
      for (i in seq_along(x$goodness.of.fit)) {
        cat("  ", goodness_names[[i]], ": ",
            format(unname(x$goodness.of.fit[[i]]), digits = 4), "\n", sep = "")
      }
    }
  }
  invisible(x)
}

.compute_wdstest_goodness <- function(dm, f, formula = NULL, formula_data = NULL,
                                      dm_adjustment_residual = NULL,
                                      goodness = "auto",
                                      distance.diagnostics = TRUE,
                                      keep.eigenvalues = FALSE, ...) {
  has_formula <- !is.null(formula)
  goodness <- .resolve_goodness(goodness, has_formula = has_formula)
  diagnostics <- list()

  if (has_formula && isTRUE(distance.diagnostics)) {
    diagnostics[["adjustment"]] <- .collect_distance_diagnostics(
      dm_adjustment_residual,
      role = "adjustment",
      formula = paste(deparse(formula), collapse = "")
    )
  }
  if (!has_formula && isTRUE(distance.diagnostics) &&
      !is.null(attr(dm, "distance.diagnostics", exact = TRUE))) {
    diagnostics[["input"]] <- .collect_distance_diagnostics(
      dm,
      role = "input"
    )
  }

  if (identical(goodness$requested, "none")) {
    return(list(
      goodness.of.fit = NULL,
      goodness.components = NULL,
      distance.diagnostics = .combine_distance_diagnostics(diagnostics)
    ))
  }

  # Treat f as an internal formula variable so users do not need to add the
  # tested factor to formula_data themselves. Name generation avoids collisions
  # with existing columns or environment bindings.
  factor_name <- .make_internal_factor_name(formula_data)
  factor_data <- .add_internal_factor(formula_data, f, factor_name)
  factor_formula <- stats::as.formula(paste("~", factor_name))
  components <- list()

  if (has_formula) {
    # The adjustment residual matrix is already produced for the adjusted test.
    r2_adjustment <- unname(dist.goodness.of.fit(dm, dm_adjustment_residual))
    components[["adjustment"]] <- .goodness_component(
      component = "adjustment",
      pseudo.R2 = r2_adjustment,
      formula = paste(deparse(formula), collapse = ""),
      interpretation = "Distance variation explained by the adjustment formula alone, ignoring the tested factor."
    )

    # factor.only requires its own F-only residual distance matrix. It is
    # skipped for adjusted goodness="auto" to avoid unnecessary computation.
    needs_factor_only <- any(goodness$computed %in% "factor.only")
    if (needs_factor_only) {
      dm_factor_residual <- .a_dist_for_goodness(
        dm,
        factor_formula,
        factor_data,
        ...,
        distance.diagnostics = distance.diagnostics,
        keep.eigenvalues = keep.eigenvalues
      )
      diagnostics[["factor.only"]] <- .collect_distance_diagnostics(
        dm_factor_residual,
        role = "factor.only",
        formula = "~ f"
      )
      r2_factor <- unname(dist.goodness.of.fit(dm, dm_factor_residual))
      components[["factor.only"]] <- .goodness_component(
        component = "factor.only",
        pseudo.R2 = r2_factor,
        formula = "~ f",
        interpretation = "Distance variation explained by the tested factor alone, ignoring adjustment variables."
      )
    }

    needs_full <- any(goodness$computed %in% c("full", "semi.partial", "partial"))
    if (needs_full) {
      # Fit Z + F from the original distances, avoiding sequential residualization.
      full_formula <- .add_factor_to_formula(formula, factor_name, formula_data)
      dm_full_residual <- .a_dist_for_goodness(
        dm,
        full_formula,
        factor_data,
        ...,
        distance.diagnostics = distance.diagnostics,
        keep.eigenvalues = keep.eigenvalues
      )
      diagnostics[["full"]] <- .collect_distance_diagnostics(
        dm_full_residual,
        role = "full",
        formula = .display_formula_with_factor(formula, formula_data)
      )
      r2_full <- unname(dist.goodness.of.fit(dm, dm_full_residual))
      r2_semi_partial <- r2_full - r2_adjustment
      adjustment_residual_fraction <- 1 - r2_adjustment

      # partial R2 is undefined if the adjustment model leaves no residual
      # variation to be explained by the tested factor.
      r2_partial <- if (is.na(r2_adjustment) ||
                        abs(adjustment_residual_fraction) < .Machine$double.eps) {
        NA_real_
      } else {
        r2_semi_partial / adjustment_residual_fraction
      }

      components[["full"]] <- .goodness_component(
        component = "full",
        pseudo.R2 = r2_full,
        formula = .display_formula_with_factor(formula, formula_data),
        interpretation = "Distance variation explained jointly by the adjustment formula and the tested factor."
      )
      components[["semi.partial"]] <- .goodness_component(
        component = "semi.partial",
        pseudo.R2 = r2_semi_partial,
        formula = "f | adjustment",
        interpretation = "Additional fraction of total distance variation explained by the tested factor after adjustment."
      )
      components[["partial"]] <- .goodness_component(
        component = "partial",
        pseudo.R2 = r2_partial,
        formula = "f | adjustment",
        interpretation = "Fraction of adjustment-residual distance variation explained by the tested factor."
      )
    }
  } else {
    if (any(goodness$computed %in% "adjustment")) {
      stop("'adjustment' goodness-of-fit requires an adjustment 'formula'.")
    }

    dm_factor_residual <- .a_dist_for_goodness(
      dm,
      factor_formula,
      factor_data,
      ...,
      distance.diagnostics = distance.diagnostics,
      keep.eigenvalues = keep.eigenvalues
    )
    diagnostics[["factor.only"]] <- .collect_distance_diagnostics(
      dm_factor_residual,
      role = "factor.only",
      formula = "~ f"
    )
    r2_factor <- unname(dist.goodness.of.fit(dm, dm_factor_residual))

    # With no adjustment set, full, semi-partial, and partial reduce to the
    # tested factor's model pseudo-R-squared.
    for (component in goodness$computed) {
      components[[component]] <- .goodness_component(
        component = component,
        pseudo.R2 = r2_factor,
        formula = "~ f",
        interpretation = .unadjusted_goodness_interpretation(component)
      )
    }
  }

  goodness.components <- do.call(rbind, components)
  goodness.components$displayed <- goodness.components$component %in% goodness$displayed
  rownames(goodness.components) <- NULL

  displayed_index <- match(goodness$displayed, goodness.components$component)
  goodness.of.fit <- goodness.components$pseudo.R2[displayed_index]
  names(goodness.of.fit) <- .goodness_component_label(goodness$displayed, has_formula)

  list(
    goodness.of.fit = goodness.of.fit,
    goodness.components = goodness.components,
    distance.diagnostics = .combine_distance_diagnostics(diagnostics)
  )
}

.resolve_goodness <- function(goodness, has_formula) {
  allowed <- c("auto", "factor.only", "adjustment", "full",
               "semi.partial", "partial", "all", "none")
  goodness <- unname(as.character(goodness))
  if (length(goodness) < 1) {
    stop("'goodness' must contain at least one value.")
  }
  if (anyNA(goodness)) {
    stop("'goodness' cannot contain missing values.")
  }
  unknown <- setdiff(goodness, allowed)
  if (length(unknown) > 0) {
    stop("Unknown 'goodness' value: ", paste(unknown, collapse = ", "), ".")
  }
  if (length(intersect(goodness, c("auto", "all", "none"))) > 0 && length(goodness) > 1) {
    stop("'auto', 'all', and 'none' cannot be combined with other 'goodness' values.")
  }

  if (identical(goodness, "none")) {
    return(list(requested = "none", displayed = "none", computed = "none"))
  }
  if (identical(goodness, "auto")) {
    if (has_formula) {
      # Adjusted auto shows the factor's semi-partial pseudo-R2. The stored
      # components include partial because it is derived from the same scalars.
      return(list(
        requested = "auto",
        displayed = "semi.partial",
        computed = c("adjustment", "full", "semi.partial", "partial")
      ))
    }
    # With no adjustment variables, the natural default is the factor-only
    # pseudo-R2 rather than a semi-partial label.
    return(list(
      requested = "auto",
      displayed = "factor.only",
      computed = "factor.only"
    ))
  }
  if (identical(goodness, "all")) {
    if (has_formula) {
      all_components <- c("adjustment", "factor.only", "full", "semi.partial", "partial")
    } else {
      all_components <- c("factor.only", "full", "semi.partial", "partial")
    }
    return(list(requested = "all", displayed = all_components, computed = all_components))
  }

  computed <- unique(goodness)
  if (has_formula) {
    # semi.partial and partial both require the adjustment-only and full-model
    # pseudo-R2 values, so store those intermediate components too.
    if (any(computed %in% c("semi.partial", "partial"))) {
      computed <- unique(c("adjustment", "full", "semi.partial", "partial", computed))
    }
    if (any(computed %in% "full")) {
      computed <- unique(c("adjustment", computed))
    }
    if (any(computed %in% "factor.only")) {
      computed <- unique(c("adjustment", computed))
    }
  }

  list(requested = goodness, displayed = goodness, computed = computed)
}

.a_dist_for_goodness <- function(dm, formula, formula_data, ...,
                                 distance.diagnostics = TRUE,
                                 keep.eigenvalues = FALSE) {
  # Keep WdS.test output focused; direct a.dist() calls still report eigenvalue
  # handling messages to users.
  suppressMessages(a.dist(
    dm,
    formula,
    formula_data,
    ...,
    distance.diagnostics = distance.diagnostics,
    keep.eigenvalues = keep.eigenvalues
  ))
}

.collect_distance_diagnostics <- function(dm, role, formula = NULL) {
  diagnostics <- attr(dm, "distance.diagnostics", exact = TRUE)
  if (is.null(diagnostics)) {
    return(NULL)
  }

  # Roles are assigned by WdS.test because a.dist() does not know whether its
  # residual matrix is being used for adjustment, factor-only, or full-model work.
  diagnostics$role <- role
  if (!is.null(formula)) {
    diagnostics$formula <- formula
  }
  diagnostics
}

.combine_distance_diagnostics <- function(diagnostics) {
  diagnostics <- Filter(Negate(is.null), diagnostics)
  if (length(diagnostics) == 0) {
    return(NULL)
  }

  combined <- do.call(rbind, diagnostics)
  rownames(combined) <- NULL
  combined
}

.make_internal_factor_name <- function(data, base_name = ".wdstar_f") {
  if (is.null(data)) {
    return(base_name)
  }

  name_exists <- function(candidate) {
    if (is.environment(data)) {
      exists(candidate, envir = data, inherits = FALSE)
    } else {
      candidate %in% names(data)
    }
  }

  candidate <- base_name
  counter <- 1L
  while (name_exists(candidate)) {
    candidate <- paste0(base_name, counter)
    counter <- counter + 1L
  }
  candidate
}

.add_internal_factor <- function(data, f, factor_name) {
  if (is.null(data)) {
    data <- data.frame(f, check.names = FALSE)
    names(data) <- factor_name
    return(data)
  }
  if (is.environment(data)) {
    factor_env <- new.env(parent = data)
    assign(factor_name, f, envir = factor_env)
    return(factor_env)
  }

  data[[factor_name]] <- f
  data
}

.add_factor_to_formula <- function(formula, factor_name, data) {
  formula_terms <- stats::terms(formula, data = data)
  term_labels <- attr(formula_terms, "term.labels")
  rhs <- paste(c(term_labels, factor_name), collapse = " + ")

  # Preserve explicit no-intercept formulas. Without this, ~ x - 1 would become
  # ~ x + f when building the full model, silently reintroducing an intercept
  # and changing full, semi.partial, and partial pseudo-R2 values.
  if (identical(attr(formula_terms, "intercept"), 0L)) {
    rhs <- paste(rhs, "- 1")
  }
  stats::as.formula(paste("~", rhs), env = environment(formula))
}

.display_formula_with_factor <- function(formula, data, factor_label = "f") {
  formula_terms <- stats::terms(formula, data = data)
  term_labels <- attr(formula_terms, "term.labels")
  rhs <- paste(c(term_labels, factor_label), collapse = " + ")

  # Match .add_factor_to_formula() so goodness.components displays the same
  # intercept convention used in the internal full-model calculation.
  if (identical(attr(formula_terms, "intercept"), 0L)) {
    rhs <- paste(rhs, "- 1")
  }
  paste("~", rhs)
}

.goodness_component <- function(component, pseudo.R2, formula, interpretation) {
  data.frame(
    component = component,
    pseudo.R2 = pseudo.R2,
    formula = formula,
    interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}

.goodness_component_label <- function(component, has_formula) {
  vapply(component, function(x) {
    switch(
      x,
      "adjustment" = "distance-based pseudo-R\u00B2 for adjustment variables",
      "factor.only" = "distance-based pseudo-R\u00B2 for factor only",
      "full" = "distance-based pseudo-R\u00B2 for full model",
      "semi.partial" = if (has_formula) {
        "semi-partial distance-based pseudo-R\u00B2 for factor given adjustment"
      } else {
        "semi-partial distance-based pseudo-R\u00B2 for factor"
      },
      "partial" = if (has_formula) {
        "partial distance-based pseudo-R\u00B2 for factor given adjustment"
      } else {
        "partial distance-based pseudo-R\u00B2 for factor"
      },
      x
    )
  }, character(1))
}

.unadjusted_goodness_interpretation <- function(component) {
  switch(
    component,
    "factor.only" = "Distance variation explained by the tested factor.",
    "full" = "Distance variation explained by the full model; this equals factor.only because no adjustment formula was supplied.",
    "semi.partial" = "Additional fraction of total distance variation explained by the tested factor; this equals factor.only because no adjustment formula was supplied.",
    "partial" = "Fraction of residual distance variation explained by the tested factor; this equals factor.only because no adjustment formula was supplied.",
    "Distance variation explained by the requested component."
  )
}

.format_estimate_name <- function(name,
                                  utf8 = isTRUE(l10n_info()[["UTF-8"]]),
                                  html = .is_html_output(),
                                  latex = .is_latex_output()) {
  if (.use_ascii_math_labels(utf8 = utf8, html = html, latex = latex)) {
    return(sub("\u03C9\u00B2", "omega^2", name, fixed = TRUE))
  }
  name
}

.format_r2_name <- function(name,
                            utf8 = isTRUE(l10n_info()[["UTF-8"]]),
                            html = .is_html_output(),
                            latex = .is_latex_output()) {
  if (.use_ascii_math_labels(utf8 = utf8, html = html, latex = latex)) {
    return(sub("R\u00B2", "R^2", name, fixed = TRUE))
  }
  name
}

.use_ascii_math_labels <- function(
    utf8 = isTRUE(l10n_info()[["UTF-8"]]),
    html = .is_html_output(),
    latex = .is_latex_output()) {
  isTRUE(latex) || (!isTRUE(utf8) && isTRUE(html))
}

.is_html_output <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    isTRUE(knitr::is_html_output())
  } else {
    FALSE
  }
}

.is_latex_output <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    isTRUE(knitr::is_latex_output())
  } else {
    FALSE
  }
}
