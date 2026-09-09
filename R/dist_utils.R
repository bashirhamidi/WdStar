#' Test if an object is of class 'dist'
#'
#' This function checks whether the given object belongs to class 'dist'.
#'
#' @param x Object to test.
#'
#' @return Logical indicating whether the object is a distance matrix.
#' @export
#' @examples
#'
#' is.dist(as.dist(matrix(1:4, nrow = 2)))
#' is.dist(matrix(1:4, nrow = 2))
#'
is.dist <- function(x) any(class(x) == "dist")

.as_formula_data <- function(formula_data) {
  if (is.environment(formula_data)) {
    return(formula_data)
  }

  data <- tryCatch(
    data.frame(formula_data, check.names = FALSE),
    error = function(e) e
  )

  if (inherits(data, "error")) {
    stop(
      "'formula_data' must be an environment, data frame, list, or coercible to a data frame.",
      call. = FALSE
    )
  }

  data
}

#' Calculate Sigma Squared for Distance Matrix
#'
#' This function computes the sigma squared statistic for a given distance matrix.
#'
#' @param dm Distance matrix.
#'
#' @return Sigma squared value.
#' @export
#' @examples
#'
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'
#'   dist.sigma2(dm)
#' }
#'
dist.sigma2 <- function(dm) {
  dd <- as.matrix(dm)
  dd[upper.tri(dd)] <- 0
  sum(dd^2) / nrow(dd) / (nrow(dd) - 1)
}

#' Calculate Distance-Based Pseudo-R-Squared for Residual Distances
#'
#' This function computes a distance-based pseudo-coefficient of determination
#' (pseudo-\eqn{R^2}) by comparing total variation in a raw distance matrix with
#' residual variation in another distance matrix.
#'
#' @param dm The original/raw distance matrix.
#' @param dm_residual A residual distance matrix after removing a model
#'   component, such as a grouping factor, adjustment variables, or their full
#'   model.
#' @param adjusted_dm Optional compatibility alias for \code{dm_residual}.
#'
#' @return Distance-based pseudo-\eqn{R^2}. Returns \code{NA_real_} when the
#'   original distance matrix has zero total variation.
#' @details Let \eqn{V(D) = [N(N-1)]^{-1}\sum_{i<j} d_{ij}^2}. The returned
#'   diagnostic is
#'   \deqn{R^2_{pseudo} = 1 - V(D_{residual}) / V(D_{original}).}
#'   The shared sample-size factor cancels in the ratio. This is a
#'   distance-based variation-reduction diagnostic, not ordinary least-squares
#'   \eqn{R^2} for a univariate response. Values near zero indicate little
#'   reduction in distance variation, larger positive values indicate greater
#'   reduction, and negative values may occur when the residual distance matrix
#'   has more variation than the original distance matrix.
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'
#'   dm_residual <- a.dist(
#'     dm = dm,
#'     formula = ~ SeqTech,
#'     formula_data = phyloseq::sample_data(ent)
#'   )
#'   dist.goodness.of.fit(dm, dm_residual)
#' }
#'
dist.goodness.of.fit <- function(dm, dm_residual = NULL, adjusted_dm = NULL) {
  # adjusted_dm is retained as a compatibility alias from the first development
  # version of this helper. New code should prefer dm_residual because the
  # residual matrix may come from a factor, adjustment variables, or a full model.
  if (!is.null(adjusted_dm)) {
    if (!is.null(dm_residual)) {
      stop("Provide only one of 'dm_residual' or 'adjusted_dm'.")
    }
    dm_residual <- adjusted_dm
  }
  if (is.null(dm_residual)) {
    stop("'dm_residual' must be provided.")
  }
  if (!is.dist(dm) || !is.dist(dm_residual)) {
    stop("'dm' and 'dm_residual' must both be distance matrices of class 'dist'.")
  }
  if (attr(dm, "Size") != attr(dm_residual, "Size")) {
    stop("'dm' and 'dm_residual' must contain the same number of observations.")
  }

  # dist.sigma2() is proportional to distance-based total SS; the shared
  # sample-size scaling cancels in SS_residual / SS_total.
  ss_total <- dist.sigma2(dm)
  ss_residual <- dist.sigma2(dm_residual)

  # If all original distances are zero, the total variation denominator is zero,
  # so the pseudo-R2 is undefined rather than 0.
  goodness.of.fit <- if (ss_total == 0) NA_real_ else 1 - (ss_residual / ss_total)
  attr(goodness.of.fit, "names") <- "distance-based pseudo-R\u00B2"
  goodness.of.fit
}

#' Calculate Sum of Squares for Distance Matrix with Factor
#'
#' This function calculates the sum of squares for a given distance matrix,
#' weighted by a factor variable.
#'
#' @param dm2 Square distance matrix.
#' @param f Factor variable for weighting.
#'
#' @return Sum of squares matrix.
#' @export
#' @examples
#'
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   dist.ss2(as.matrix(dm)^2, f)
#' }
#'
dist.ss2 <- function(dm2, f) {
  K <- sapply(levels(f), function(lev) f == lev)
  t(K) %*% dm2 %*% K / 2
}

#' Calculate Group-wise Sigma Squared
#'
#' This function calculates the sigma squared statistic for each group defined
#' by a factor variable in a given distance matrix.
#'
#' @param dm Distance matrix.
#' @param f Factor variable for group definition.
#'
#' @return A named numeric vector of group-wise sigma squared values.
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   dist.group.sigma2(dm, f)
#' }
dist.group.sigma2 <- function(dm, f) {
  diag(dist.ss2(as.matrix(dm)^2, f)) / table(f) / (table(f) - 1)
}

#' Calculate Cohen's d for Distance Matrix
#'
#' This function calculates the Cohen's d statistic for a distance matrix,
#' assuming two levels in the factor variable.
#'
#' @param dm Distance matrix.
#' @param f Factor variable with two levels.
#'
#' @return Cohen's d value if factor has exactly two levels; NULL otherwise.
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   keep <- f %in% c("1", "2")
#'   dm_12 <- as.dist(as.matrix(dm)[keep, keep])
#'   f_12 <- droplevels(f[keep])
#'
#'   dist.cohen.d(dm_12, f_12)
#' }
#'
dist.cohen.d <- function(dm, f) {
  if (nlevels(f) != 2) {
    return(NULL)
  }
  SS2 <- dist.ss2(as.matrix(dm^2), f)
  ns <- summary(f)
  N <- sum(ns)

  SST <- sum(SS2) / N
  SSW <- sum(diag(SS2) / ns)

  mean.diff <- (sqrt((ns[1] + ns[2]) / (ns[1] * ns[2]))) * sqrt(SST - SSW)

  sigmas <- diag(SS2) / ns / (ns - 1)
  s1 <- sigmas[1]
  s2 <- sigmas[2]

  mean.diff / sqrt(((ns[1] - 1) * s1 + (ns[2] - 1) * s2) / (sum(ns) - 2))
}

#' Compute a Covariate-Adjusted Distance Matrix
#'
#' This function takes a formula and distance matrix, projects out covariate
#' effects through principal-coordinate residualization, and returns the
#' resulting adjusted distance matrix.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#'
#' @param formula Only the right-hand side of a typical formula, such as
#'   \code{~ SeqTech} or \code{~ Age + SeqTech}, is necessary. The formula has
#'   similar requirements as in \code{vegan::adonis()}.
#'
#' @param formula_data A dataset which contains the variables specified in
#'   formula. It may be an environment, data frame, list, or object coercible to
#'   a data frame, such as \code{phyloseq::sample_data()}. Row names should match
#'   the row names in distance matrix \code{dm}. If not provided, the parent
#'   frame will be used.
#' @param tol Tolerance for eigenvalues. This is the cutoff for the eigenvalues
#'   to be considered zero. Default is 10^-8.
#' @param distance.diagnostics Logical indicating whether to attach a compact
#'   diagnostic summary of the eigendecomposition used to reconstruct the
#'   residual distance matrix. Default is \code{TRUE}.
#' @param keep.eigenvalues Logical indicating whether to store the raw
#'   eigenvalues in the diagnostic summary. Default is \code{FALSE} to keep
#'   repeated simulations and large analyses memory-light.
#'
#' @return A covariate-adjusted distance matrix of class \code{dist}. When
#'   \code{distance.diagnostics = TRUE}, the returned object has a
#'   \code{"distance.diagnostics"} attribute containing one row of eigenvalue
#'   diagnostics.
#'
#' @details The \code{a.dist()} function only requires a right-hand side of the
#'   formula. Instead of the left-hand side, it uses the dissimilarity distance
#'   matrix \code{dm}. The function constructs a model matrix from the
#'   right-hand side (RHS) of the formula. After performing necessary matrix
#'   operations and eigen-decomposition, it calculates the Euclidean distances.
#'   It preserves the labels of the input dm.
#'
#'   During reconstruction, eigenvalues with \code{abs(lambda) < tol} are
#'   treated as zero, and remaining negative eigenvalues are set to zero before
#'   distances are recomputed. The \code{"distance.diagnostics"} attribute records
#'   how many eigenvalues were positive, negative, treated as zero by
#'   \code{tol}, and removed as negative after tolerance handling. These
#'   diagnostics describe construction of the residual distance matrix and are
#'   separate from goodness-of-fit statistics such as pseudo-\eqn{R^2}. If
#'   \code{keep.eigenvalues = TRUE}, the diagnostics also include the raw
#'   eigenvalues before tolerance and negative-value zeroing.
#'
#'   The diagnostic fields include:
#'   \itemize{
#'     \item \code{role}: how \code{WdS.test()} used the matrix, if applicable;
#'     direct \code{a.dist()} calls leave this as \code{NA}.
#'     \item \code{formula} and \code{tol}: the formula and tolerance used.
#'     \item \code{n.eigenvalues}: number of eigenvalues from the decomposition.
#'     \item \code{n.positive.raw}, \code{n.negative.raw}, and
#'     \code{n.zero.raw}: signs before tolerance handling.
#'     \item \code{n.zeroed.by.tol}: eigenvalues with
#'     \code{abs(lambda) < tol}, including exact zeros.
#'     \item \code{n.positive.zeroed.by.tol} and
#'     \code{n.negative.zeroed.by.tol}: positive and negative raw eigenvalues
#'     counted in \code{n.zeroed.by.tol}.
#'     \item \code{n.negative.removed}: negative eigenvalues remaining after
#'     tolerance handling, which are set to zero before reconstruction.
#'     \item \code{n.positive.retained} and \code{n.zero.final}: counts after
#'     tolerance and negative-value zeroing.
#'     \item \code{min.eigenvalue}, \code{max.eigenvalue},
#'     \code{sum.positive.raw}, \code{sum.negative.raw},
#'     \code{sum.negative.magnitude.raw}, and \code{prop.negative.raw}: compact
#'     summaries of the raw eigenspectrum.
#'     \item \code{eigenvalues}: optional raw eigenvalue vector, present only
#'     when \code{keep.eigenvalues = TRUE}.
#'   }
#'
#'   This function refactors and generalizes functionality from
#'   \code{aPCoA::aPCoA()} function in the aPCoA package.
#' @references  Shi Y, Zhang L, Do KA, Peterson CB, Jenq RR. aPCoA: covariate
#'   adjusted principal coordinates analysis. Bioinformatics.
#'   2020;36(13):4099-4101. doi:10.1093/bioinformatics/btaa276
#'
#'   Shi Y (2021). aPCoA: Covariate Adjusted PCoA Plot. R package version 1.3,
#'   https://CRAN.R-project.org/package=aPCoA.
#'
#'   Please cite both the package and the paper when using this function.
#'
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'
#'   ## Use the full enterotype data after removing samples without Enterotype.
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'
#'   ## Adjust the distance matrix by sequencing technology, a categorical
#'   ## covariate stored in phyloseq::sample_data(ent).
#'   a.dm <- a.dist(
#'     dm = dm,
#'     formula = ~ SeqTech,
#'     formula_data = phyloseq::sample_data(ent)
#'   )
#'   a.dm
#'   attr(a.dm, "distance.diagnostics")
#'
#'   ## Continuous covariates can also be used. Age has many missing values in
#'   ## enterotype, so this example builds a matching Age-complete distance
#'   ## matrix before using Age in the adjustment formula.
#'   ent_age <- phyloseq::subset_samples(
#'     enterotype,
#'     !is.na(Enterotype) & !is.na(Age)
#'   )
#'   dm_age <- phyloseq::distance(ent_age, method = "bray")
#'   a.dm.age <- a.dist(
#'     dm = dm_age,
#'     formula = ~ Age,
#'     formula_data = phyloseq::sample_data(ent_age)
#'   )
#'   attr(a.dm.age, "distance.diagnostics")
#' }
#'
a.dist = function(dm, formula, formula_data=parent.frame(), tol=10^-8,
                  distance.diagnostics = TRUE, keep.eigenvalues = FALSE)
{
  .validate_distance_diagnostic_options(distance.diagnostics, keep.eigenvalues)
  if (!is.numeric(tol) || length(tol) != 1 || is.na(tol) || tol < 0) {
    stop("'tol' must be a single non-negative numeric value.")
  }

  data <- .as_formula_data(formula_data)
  Terms <- stats::terms(formula, data = data)
  # lhs <- formula[[2]]
  # lhs <- eval(lhs, data, parent.frame())
  # formula[[2]] <- NULL
  lhs <- dm
  rhs.frame <- stats::model.frame(formula, data, drop.unused.levels = TRUE) ## extract data.frame from a formula
  rhs <- stats::model.matrix(formula, rhs.frame) ## expand and create dummy variables for the terms
  grps <- attr(rhs, "assign") ## tells which column from the data set each dummy variable came from
  qrhs <- qr(rhs) ##computes the QR decomposition of rhs matrix
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1)
    stop("right-hand-side of formula has no usable terms")
  dmat <- as.matrix(lhs^2) #still has the labels of lhs
  X <- rhs
  X <- as.matrix(X)
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  A <- -1/2 * dmat
  J <- diag(nrow(X)) - matrix(rep(1/(nrow(X)), length(A)), nrow = nrow(A))
  E <- (diag(nrow(H)) - H) %*% J %*% A %*% J %*% (diag(nrow(H)) - H)

  # rownames(E) <- rownames(data)
  # colnames(E) <- rownames(data)

  # Correct rounding errors that cause E to be not symmetric
  E <- (E + t(E))/2

  eig <- eigen(E)
  lambda_raw <- eig$values
  lambda <- lambda_raw

  below_tol <- abs(lambda_raw) < tol
  if (any(below_tol)) {
    message(paste(
      sum(below_tol), "out of", length(lambda),
      "eigenvalues were smaller than 'tol' and were set to zero."
    ))
    if (missing(tol)) {
      message("The default 'tol' is 1e-8; use a smaller 'tol' to retain more small eigenvalues.")
    }
    lambda[below_tol] <- 0
  }

  negative <- lambda < 0
  if (any(negative)) {
    message(paste(
      sum(negative), "out of", length(lambda),
      "eigenvalues were negative and were set to zero."
    ))
    lambda[negative] <- 0
  }

  w <- t(t(eig$vectors) * sqrt(lambda))
  w <- stats::dist(w)
  attr(w, "Labels") <- attr(lhs, "Labels")
  if (isTRUE(distance.diagnostics)) {
    # Diagnostics are stored as metadata on the residual distance object. This
    # keeps a.dist() users informed without changing the dist values themselves.
    attr(w, "distance.diagnostics") <- .make_distance_diagnostics(
      lambda.raw = lambda_raw,
      lambda.final = lambda,
      formula = formula,
      tol = tol,
      keep.eigenvalues = keep.eigenvalues
    )
  }
  return(w)
}

.validate_distance_diagnostic_options <- function(distance.diagnostics, keep.eigenvalues) {
  if (!is.logical(distance.diagnostics) ||
      length(distance.diagnostics) != 1 ||
      is.na(distance.diagnostics)) {
    stop("'distance.diagnostics' must be TRUE or FALSE.")
  }
  if (!is.logical(keep.eigenvalues) ||
      length(keep.eigenvalues) != 1 ||
      is.na(keep.eigenvalues)) {
    stop("'keep.eigenvalues' must be TRUE or FALSE.")
  }
  if (isTRUE(keep.eigenvalues) && !isTRUE(distance.diagnostics)) {
    stop("'keep.eigenvalues = TRUE' requires 'distance.diagnostics = TRUE'.")
  }
}

.make_distance_diagnostics <- function(lambda.raw, lambda.final, formula, tol,
                                       keep.eigenvalues = FALSE) {
  raw.positive <- lambda.raw > 0
  raw.negative <- lambda.raw < 0
  raw.zero <- lambda.raw == 0
  below.tol <- abs(lambda.raw) < tol
  total.magnitude <- sum(abs(lambda.raw))
  negative.magnitude <- sum(abs(lambda.raw[raw.negative]))

  diagnostics <- data.frame(
    role = NA_character_,
    formula = paste(deparse(formula), collapse = ""),
    tol = tol,
    n.eigenvalues = length(lambda.raw),
    n.positive.raw = sum(raw.positive),
    n.negative.raw = sum(raw.negative),
    n.zero.raw = sum(raw.zero),
    n.zeroed.by.tol = sum(below.tol),
    n.positive.zeroed.by.tol = sum(raw.positive & below.tol),
    n.negative.zeroed.by.tol = sum(raw.negative & below.tol),
    n.negative.removed = sum(raw.negative & !below.tol),
    n.positive.retained = sum(lambda.final > 0),
    n.zero.final = sum(lambda.final == 0),
    min.eigenvalue = min(lambda.raw),
    max.eigenvalue = max(lambda.raw),
    sum.positive.raw = sum(lambda.raw[raw.positive]),
    sum.negative.raw = sum(lambda.raw[raw.negative]),
    sum.negative.magnitude.raw = negative.magnitude,
    prop.negative.raw = if (total.magnitude == 0) NA_real_ else negative.magnitude / total.magnitude,
    keep.eigenvalues = isTRUE(keep.eigenvalues),
    stringsAsFactors = FALSE
  )

  if (isTRUE(keep.eigenvalues)) {
    # Store raw eigenvalues only when requested. They can be useful for deep
    # diagnosis but add memory cost when many tests are retained.
    diagnostics$eigenvalues <- I(list(lambda.raw))
  }

  diagnostics
}
