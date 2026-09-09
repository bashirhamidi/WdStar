#' Compute Leave-One-Feature-Out Importance Scores
#'
#' This function estimates the importance of each adjustment feature in a
#' \code{\link{WdS.test}} model by removing one term from the adjustment formula
#' at a time. For each removed term, it reruns \code{\link{WdS.test}} with the
#' remaining terms and ranks features by the loss in adjustment goodness-of-fit.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param formula A right-hand side ONLY formula containing the adjustment
#'   features to evaluate.
#' @param formula_data An environment, data frame, list, or object coercible to a
#'   data frame, such as \code{phyloseq::sample_data()}, to be used with
#'   \code{formula}. Default is \code{parent.frame()}.
#' @param nrep The number of permutations to conduct for each
#'   \code{\link{WdS.test}} call. Default is 999.
#' @param strata A factor variable representing strata. Default is \code{NULL}.
#' @param sort Logical indicating whether to sort rows by importance. Default is
#'   \code{TRUE}.
#' @param decreasing Logical indicating the sort direction when \code{sort} is
#'   \code{TRUE}. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{WdS.test}}, such as
#'   \code{tol} for eigenvalue handling in \code{\link{a.dist}}.
#'
#' @return A data frame with one row per adjustment feature. The
#'   \code{importance} column is the decrease in adjustment-only
#'   distance-based pseudo-\eqn{R^2} after removing the feature. The
#'   \code{test.result} column
#'   stores the full \code{\link{WdS.test}} result from the leave-one-feature-out
#'   model.
#'
#' @details Importance is computed as the difference between adjustment-only
#' pseudo-\eqn{R^2} for the complete adjustment formula and adjustment-only
#' pseudo-\eqn{R^2} after removing one feature:
#' \deqn{R^2_{pseudo,all\ adjustment\ features} -
#'       R^2_{pseudo,without\ feature}.}
#' Larger values indicate that removing the feature caused a larger decrease in
#' the variation explained by the adjustment formula.
#'
#' @seealso \code{\link{WdS.test}}, \code{\link{dist.goodness.of.fit}}
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'
#'   ## Age is missing for many enterotype samples, so use a matching
#'   ## complete-case distance matrix when ranking Age and Gender.
#'   ent <- phyloseq::subset_samples(
#'     enterotype,
#'     !is.na(Enterotype) & !is.na(Age) & !is.na(Gender)
#'   )
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   WdS.feature.importance(
#'     dm = dm,
#'     f = f,
#'     formula = ~ Gender + Age,
#'     formula_data = phyloseq::sample_data(ent),
#'     nrep = 9
#'   )
#' }
WdS.feature.importance <- function(dm, f, formula, formula_data = parent.frame(),
                                   nrep = 999, strata = NULL, sort = TRUE,
                                   decreasing = TRUE, ...) {
  if (missing(formula) || is.null(formula)) {
    stop("'formula' must be provided to compute feature importance.")
  }

  data <- .as_formula_data(formula_data)
  formula_terms <- attr(stats::terms(formula, data = data), "term.labels")

  if (length(formula_terms) < 1) {
    stop("'formula' must contain at least one adjustment feature.")
  }

  full_result <- WdS.test(
    dm = dm,
    f = f,
    nrep = nrep,
    strata = strata,
    formula = formula,
    formula_data = data,
    # Feature importance ranks adjustment features, so use the adjustment-only
    # pseudo-R2 rather than WdS.test()'s default factor/semi-partial component.
    goodness = "adjustment",
    ...
  )
  full_goodness_of_fit <- unname(full_result$goodness.of.fit)

  rows <- vector("list", length(formula_terms))
  test_results <- vector("list", length(formula_terms))

  for (i in seq_along(formula_terms)) {
    removed_feature <- formula_terms[[i]]
    remaining_features <- setdiff(formula_terms, removed_feature)

    # With one original formula term, removing it leaves the unadjusted WdS.test.
    if (length(remaining_features) == 0) {
      leave_one_out_formula <- NULL
      leave_one_out_formula_text <- NA_character_
      leave_one_out_result <- WdS.test(
        dm = dm,
        f = f,
        nrep = nrep,
        strata = strata,
        # Removing the only adjustment feature leaves no adjustment model; use
        # zero as the adjustment-only goodness baseline.
        goodness = "none"
      )
      leave_one_out_goodness_of_fit <- 0
    } else {
      leave_one_out_formula <- stats::as.formula(
        paste("~", paste(remaining_features, collapse = " + ")),
        env = environment(formula)
      )
      leave_one_out_formula_text <- paste(deparse(leave_one_out_formula), collapse = "")
      leave_one_out_result <- WdS.test(
        dm = dm,
        f = f,
        nrep = nrep,
        strata = strata,
        formula = leave_one_out_formula,
        formula_data = data,
        # Keep the leave-one-out score on the same adjustment-only scale as the
        # complete adjustment formula.
        goodness = "adjustment",
        ...
      )
      leave_one_out_goodness_of_fit <- unname(leave_one_out_result$goodness.of.fit)
    }

    importance <- full_goodness_of_fit - leave_one_out_goodness_of_fit

    rows[[i]] <- data.frame(
      feature = removed_feature,
      importance = importance,
      full.goodness.of.fit = full_goodness_of_fit,
      leave.one.out.goodness.of.fit = leave_one_out_goodness_of_fit,
      statistic = unname(leave_one_out_result$statistic),
      p.value = leave_one_out_result$p.value,
      omega.squared = unname(leave_one_out_result$estimate),
      formula = leave_one_out_formula_text,
      stringsAsFactors = FALSE
    )
    test_results[[i]] <- leave_one_out_result
  }

  importance_table <- do.call(rbind, rows)
  importance_table$test.result <- I(test_results)

  if (isTRUE(sort)) {
    importance_table <- importance_table[
      order(importance_table$importance, decreasing = decreasing, na.last = TRUE),
      ,
      drop = FALSE
    ]
  }

  rownames(importance_table) <- NULL
  importance_table
}

#' Compute One-Taxon Adjustment Importance Scores
#'
#' This function ranks ASVs, OTUs, or other taxa by adjusting for one taxon at a
#' time in \code{\link{WdS.test}}. For each taxon, the taxon abundance vector is
#' used as the adjustment formula, and the tested factor's goodness-of-fit after
#' that adjustment is compared with the unadjusted factor goodness-of-fit.
#'
#' @param dm A distance matrix (any arbitrary distance or dissimilarity metric).
#' @param f A factor variable representing the groups.
#' @param physeq Optional \code{phyloseq} object. When supplied,
#'   \code{taxa_table} is extracted with \code{phyloseq::otu_table()} and
#'   taxonomy is extracted with \code{phyloseq::tax_table()}, if available.
#' @param taxa_table Optional numeric ASV/OTU/taxa abundance table. Supply this
#'   when \code{physeq} is not provided.
#' @param taxonomy_table Optional taxa-by-rank taxonomy table. Row names should
#'   match taxon names in \code{taxa_table}. If row names are absent, rows are
#'   assumed to be in the same taxon order as \code{taxa_table}.
#' @param taxa_are_rows Logical indicating whether rows of \code{taxa_table} are
#'   taxa and columns are samples. Ignored when \code{physeq} is supplied because
#'   the orientation is read from the \code{phyloseq} object.
#' @param taxonomy.ranks Taxonomy ranks to include in the output. The default,
#'   \code{"auto"}, includes the lowest two taxonomy ranks with at least one
#'   non-missing value. Use \code{"all"} for every available rank, a character
#'   vector such as \code{c("Genus", "Species")} for selected ranks, or
#'   \code{NULL} to omit taxonomy columns.
#' @param formula Optional right-hand side formula controlling the adjustment
#'   model used for each taxon. The default \code{NULL} performs a taxon-only
#'   adjustment, equivalent to \code{~ .wdstar_taxon}. If \code{formula} does
#'   not include \code{.wdstar_taxon}, the taxon term is added automatically.
#'   For example, \code{formula = ~ Subject_ID} evaluates each taxon with
#'   \code{~ Subject_ID + .wdstar_taxon}. If \code{formula} includes
#'   \code{.wdstar_taxon}, it is used as supplied, allowing formulas such as
#'   \code{~ .wdstar_taxon * Subject_ID}.
#' @param formula_data Optional sample metadata used with \code{formula}. This
#'   may be an environment, data frame, list, or object coercible to a data
#'   frame, such as \code{phyloseq::sample_data()}. Row names are used to align
#'   data-frame-like inputs to \code{dm} labels when available.
#' @param rank.by Column used for ranking. The default,
#'   \code{"delta.factor.goodness.of.fit"}, ranks taxa by how much the tested
#'   factor's pseudo-\eqn{R^2} decreases after adjusting for that taxon.
#' @param nrep The number of permutations to conduct for each
#'   \code{\link{WdS.test}} call. Default is 999.
#' @param strata A factor variable representing strata. Default is \code{NULL}.
#' @param distance.diagnostics Logical passed to \code{\link{WdS.test}}. Default
#'   is \code{FALSE} to keep large taxa scans compact.
#' @param keep.eigenvalues Logical passed to \code{\link{WdS.test}}. Default is
#'   \code{FALSE}.
#' @param keep.tests Logical indicating whether to retain each adjusted
#'   \code{\link{WdS.test}} result in a \code{test.result} list column. Default
#'   is \code{FALSE} because this can become large for many taxa.
#' @param sort Logical indicating whether to sort rows by \code{rank.by}.
#'   Default is \code{TRUE}.
#' @param decreasing Optional logical sort direction. By default, larger values
#'   are ranked first, except for \code{rank.by = "p.value"} and
#'   \code{rank.by = "factor.goodness.of.fit.after.adjustment"}, where smaller
#'   values are ranked first.
#' @param verbose Logical indicating whether to show messages from the repeated
#'   internal \code{\link{WdS.test}} calls. Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{WdS.test}}, such as
#'   \code{tol} for eigenvalue handling in \code{\link{a.dist}}.
#'
#' @return A \code{data.table} with one row per taxon. The \code{taxon} column
#'   stores the ASV/OTU/taxon identifier, followed by the requested taxonomy
#'   ranks when available. The \code{adjustment.formula} column records the
#'   taxon-specific formula passed to \code{\link{WdS.test}}. The
#'   \code{adjustment.goodness.of.fit} column is the distance-based
#'   pseudo-\eqn{R^2} explained by that full adjustment formula. This is taxon
#'   only when \code{formula = NULL}, and taxon plus the requested sample-level
#'   terms when \code{formula} is supplied. The
#'   \code{factor.goodness.of.fit.after.adjustment} column is the tested
#'   factor's semi-partial pseudo-\eqn{R^2} after adjusting for that formula.
#'   The default \code{importance} value is
#'   \code{factor.goodness.of.fit.unadjusted -
#'   factor.goodness.of.fit.after.adjustment}.
#'
#' @details No abundance transformation is performed inside this function. Taxon
#' abundances are used exactly as supplied; transform, normalize, or aggregate
#' upstream if that is desired.
#'
#' The default ranking asks: how much of the tested factor's distance-based
#' goodness-of-fit is reduced when this taxon is adjusted out? Larger
#' \code{delta.factor.goodness.of.fit} values indicate taxa that account for more
#' of the distance structure otherwise attributed to the tested factor. The
#' separate \code{adjustment.goodness.of.fit} column answers a different
#' question: how much total distance variation is explained by the whole
#' taxon-specific adjustment formula?
#'
#' Use \code{formula} when each taxon should be evaluated with additional
#' sample-level adjustment terms. For paired or repeated-measures designs,
#' \code{formula = ~ Subject_ID} evaluates each taxon after adding
#' \code{Subject_ID} to the adjustment model. To control the exact placement of
#' the taxon term, include the reserved placeholder \code{.wdstar_taxon}
#' directly in \code{formula}.
#'
#' Constant taxa are retained in the output with \code{tested = FALSE} and
#' missing test statistics because they cannot form a usable one-taxon
#' adjustment model.
#'
#' @seealso \code{\link{WdS.test}}, \code{\link{WdS.feature.importance}},
#'   \code{\link{dist.goodness.of.fit}}
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'
#'   ## Use the full enterotype data after removing samples without Enterotype.
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   ## Rank taxa by how much the tested factor's goodness-of-fit changes after
#'   ## each taxon is adjusted out.
#'   taxa_rank <- WdS.taxa.importance(
#'     dm = dm,
#'     f = f,
#'     physeq = ent,
#'     nrep = 9
#'   )
#'   taxa_rank[, c("taxon", "Genus", "rank", "importance", "p.value"),
#'             with = FALSE]
#'
#'   ## Add a sample-level covariate to every taxon-specific adjustment model.
#'   taxa_rank_seqtech <- WdS.taxa.importance(
#'     dm = dm,
#'     f = f,
#'     physeq = ent,
#'     formula = ~ SeqTech,
#'     formula_data = phyloseq::sample_data(ent),
#'     rank.by = "adjustment.goodness.of.fit",
#'     nrep = 9
#'   )
#'   taxa_rank_seqtech[, c("taxon", "Genus", "rank", "importance"),
#'                     with = FALSE]
#' }
#' }
WdS.taxa.importance <- function(dm, f, physeq = NULL, taxa_table = NULL,
                                taxonomy_table = NULL, taxa_are_rows = TRUE,
                                taxonomy.ranks = "auto",
                                rank.by = c("delta.factor.goodness.of.fit",
                                            "adjustment.goodness.of.fit",
                                            "factor.goodness.of.fit.after.adjustment",
                                            "full.goodness.of.fit",
                                            "semi.partial.goodness.of.fit",
                                            "partial.goodness.of.fit",
                                            "delta.omega.squared",
                                            "omega.squared",
                                            "p.value",
                                            "statistic"),
                                nrep = 999, strata = NULL,
                                distance.diagnostics = FALSE,
                                keep.eigenvalues = FALSE,
                                keep.tests = FALSE,
                                sort = TRUE, decreasing = NULL,
                                verbose = FALSE, formula = NULL,
                                formula_data = parent.frame(), ...) {
  if (!is.dist(dm)) {
    stop("'dm' must be a distance matrix of class 'dist'.")
  }
  if (length(f) != attr(dm, "Size")) {
    stop("'f' must have the same length as the number of observations in 'dm'.")
  }
  .validate_logical_scalar(taxa_are_rows, "taxa_are_rows")
  .validate_logical_scalar(distance.diagnostics, "distance.diagnostics")
  .validate_logical_scalar(keep.eigenvalues, "keep.eigenvalues")
  .validate_logical_scalar(keep.tests, "keep.tests")
  .validate_logical_scalar(sort, "sort")
  .validate_logical_scalar(verbose, "verbose")
  if (!is.null(decreasing)) {
    .validate_logical_scalar(decreasing, "decreasing")
  }

  rank.by <- match.arg(rank.by)
  if (is.null(decreasing)) {
    # Most importance scores are useful with larger values first. A smaller
    # p-value and a smaller post-adjustment factor R2 are the two exceptions.
    decreasing <- !(rank.by %in% c("p.value", "factor.goodness.of.fit.after.adjustment"))
  }

  dots <- list(...)
  # WdS.taxa.importance() must control goodness itself so every taxon row has
  # comparable adjustment, semi-partial, full, and partial goodness columns.
  reserved_dots <- intersect(names(dots), "goodness")
  if (length(reserved_dots) > 0) {
    stop(
      "The following WdS.test() arguments are controlled by WdS.taxa.importance(): ",
      paste(reserved_dots, collapse = ", "),
      "."
    )
  }
  .validate_taxa_importance_formula(formula)

  # Convert either a phyloseq object or user-supplied abundance table into a
  # sample-by-taxon numeric matrix. No abundance transformation is done here.
  taxa_inputs <- .extract_taxa_importance_inputs(
    physeq = physeq,
    taxa_table = taxa_table,
    taxonomy_table = taxonomy_table,
    taxa_are_rows = taxa_are_rows
  )
  taxa_matrix <- .align_taxa_table_to_dist(taxa_inputs$taxa_table, dm)
  taxa_matrix <- .name_taxa_columns(taxa_matrix)
  if (ncol(taxa_matrix) < 1) {
    stop("'taxa_table' must contain at least one taxon.")
  }
  taxonomy_output <- .select_taxonomy_output(
    taxa_inputs$taxonomy_table,
    taxa_ids = colnames(taxa_matrix),
    taxonomy.ranks = taxonomy.ranks
  )
  taxon_summary <- .summarize_taxa_matrix(taxa_matrix)
  # formula_data is only needed when the user supplies sample-level adjustment
  # terms, such as Subject_ID. It is aligned to dm labels before the scan.
  taxa_formula_data <- .prepare_taxa_importance_formula_data(
    formula = formula,
    formula_data = formula_data,
    dm = dm
  )
  # Use one reserved placeholder for the taxon abundance vector. Real taxon IDs
  # may be numeric, contain spaces, or otherwise be invalid formula names.
  taxon_variable <- ".wdstar_taxon"
  taxon_formula <- .build_taxa_importance_formula(
    formula = formula,
    formula_data = taxa_formula_data,
    taxon_variable = taxon_variable
  )
  taxon_formula_text <- paste(deparse(taxon_formula), collapse = "")

  # The unadjusted factor-only fit does not depend on the taxon being scanned,
  # so compute it once and reuse it for delta.factor.goodness.of.fit.
  baseline_result <- .run_taxa_importance_wdstest(
    args = list(
      dm = dm,
      f = f,
      nrep = nrep,
      strata = strata,
      goodness = "factor.only",
      distance.diagnostics = distance.diagnostics,
      keep.eigenvalues = keep.eigenvalues
    ),
    dots = dots,
    verbose = verbose
  )
  factor_goodness_unadjusted <- .extract_goodness_value(baseline_result, "factor.only")
  omega_unadjusted <- unname(baseline_result$estimate)

  rows <- vector("list", ncol(taxa_matrix))
  test_results <- vector("list", ncol(taxa_matrix))

  for (i in seq_len(ncol(taxa_matrix))) {
    taxon_id <- colnames(taxa_matrix)[[i]]
    taxon_values <- taxa_matrix[, i]
    taxon_row <- taxon_summary[i, , drop = FALSE]

    # a.dist() cannot fit an adjustment term with no variation. Keep the taxon
    # in the output so users can see that it was skipped instead of silently
    # shortening or reindexing the result table.
    if (!isTRUE(taxon_row$taxon.variance > 0)) {
      rows[[i]] <- .taxa_importance_row(
        taxon = taxon_id,
        tested = FALSE,
        skip.reason = "zero variance",
        adjustment.formula = taxon_formula_text,
        factor.goodness.of.fit.unadjusted = factor_goodness_unadjusted,
        omega.squared.unadjusted = omega_unadjusted,
        taxon_summary = taxon_row
      )
      next
    }

    # The reserved .wdstar_taxon variable is replaced on each iteration. If the
    # user supplied additional sample metadata, it is carried alongside the
    # current taxon abundance vector before calling WdS.test().
    formula_data <- .add_taxon_to_formula_data(
      data = taxa_formula_data,
      taxon_values = taxon_values,
      taxon_variable = taxon_variable,
      sample_names = rownames(taxa_matrix)
    )

    # Request adjustment and semi.partial. WdS.test() will also store full and
    # partial because they are computed from the same adjustment/full-model fits.
    adjusted_result <- .run_taxa_importance_wdstest(
      args = list(
        dm = dm,
        f = f,
        nrep = nrep,
        strata = strata,
        formula = taxon_formula,
        formula_data = formula_data,
        goodness = c("adjustment", "semi.partial"),
        distance.diagnostics = distance.diagnostics,
        keep.eigenvalues = keep.eigenvalues
      ),
      dots = dots,
      verbose = verbose
    )

    adjustment_r2 <- .extract_goodness_value(adjusted_result, "adjustment")
    semi_partial_r2 <- .extract_goodness_value(adjusted_result, "semi.partial")
    full_r2 <- .extract_goodness_value(adjusted_result, "full")
    partial_r2 <- .extract_goodness_value(adjusted_result, "partial")
    omega_adjusted <- unname(adjusted_result$estimate)

    # Store several candidate ranking quantities. The importance column is just
    # an alias of rank.by, so users can re-rank this table without rerunning the
    # expensive WdS.test() loop.
    rows[[i]] <- .taxa_importance_row(
      taxon = taxon_id,
      tested = TRUE,
      skip.reason = NA_character_,
      adjustment.formula = taxon_formula_text,
      adjustment.goodness.of.fit = adjustment_r2,
      factor.goodness.of.fit.unadjusted = factor_goodness_unadjusted,
      factor.goodness.of.fit.after.adjustment = semi_partial_r2,
      delta.factor.goodness.of.fit = factor_goodness_unadjusted - semi_partial_r2,
      full.goodness.of.fit = full_r2,
      semi.partial.goodness.of.fit = semi_partial_r2,
      partial.goodness.of.fit = partial_r2,
      statistic = unname(adjusted_result$statistic),
      p.value = adjusted_result$p.value,
      omega.squared = omega_adjusted,
      omega.squared.unadjusted = omega_unadjusted,
      delta.omega.squared = omega_unadjusted - omega_adjusted,
      taxon_summary = taxon_row
    )
    test_results[[i]] <- adjusted_result
  }

  result <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (ncol(taxonomy_output) > 0) {
    for (rank_name in names(taxonomy_output)) {
      result[[rank_name]] <- taxonomy_output[[rank_name]]
    }
  }
  if (isTRUE(keep.tests)) {
    result[["test.result"]] <- test_results
  }

  # Keep a single importance column for convenient plotting/reporting while
  # retaining the original named score columns for interpretation.
  result[["importance"]] <- result[[rank.by]]
  result[["rank.by"]] <- rank.by

  if (isTRUE(sort)) {
    row_order <- order(result[[rank.by]], decreasing = decreasing, na.last = TRUE)
    # Keep row ordering explicit; data.table's single-argument `[` has special
    # evaluation rules that are unnecessary here.
    result <- result[row_order, , with = FALSE]
  }
  result[["rank"]] <- seq_len(nrow(result))

  front_columns <- c(
    "taxon",
    names(taxonomy_output),
    "rank",
    "rank.by",
    "importance",
    "tested",
    "skip.reason",
    "adjustment.formula",
    "adjustment.goodness.of.fit",
    "factor.goodness.of.fit.unadjusted",
    "factor.goodness.of.fit.after.adjustment",
    "delta.factor.goodness.of.fit",
    "full.goodness.of.fit",
    "semi.partial.goodness.of.fit",
    "partial.goodness.of.fit",
    "statistic",
    "p.value",
    "omega.squared",
    "omega.squared.unadjusted",
    "delta.omega.squared"
  )
  front_columns <- front_columns[front_columns %in% names(result)]
  final_columns <- c(front_columns, setdiff(names(result), front_columns))
  # Use data.table subsetting instead of setcolorder() so the returned object is
  # a fresh table and cannot be affected by by-reference column-order side effects.
  result <- result[, final_columns, with = FALSE]
  result
}

.extract_taxa_importance_inputs <- function(physeq, taxa_table, taxonomy_table,
                                            taxa_are_rows) {
  if (!is.null(physeq)) {
    if (!is.null(taxa_table)) {
      stop("Provide only one of 'physeq' or 'taxa_table'.")
    }
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("'physeq' input requires the phyloseq package to be installed.")
    }

    taxa_table <- as.matrix(phyloseq::otu_table(physeq))
    if (isTRUE(phyloseq::taxa_are_rows(physeq))) {
      taxa_table <- t(taxa_table)
    }

    if (is.null(taxonomy_table)) {
      physeq_taxonomy <- tryCatch(
        phyloseq::tax_table(physeq),
        error = function(e) NULL
      )
      if (!is.null(physeq_taxonomy)) {
        taxonomy_table <- as.matrix(physeq_taxonomy)
      }
    }
  } else {
    if (is.null(taxa_table)) {
      stop("Provide either 'physeq' or 'taxa_table'.")
    }
    if (isTRUE(taxa_are_rows)) {
      taxa_table <- t(as.matrix(taxa_table))
    } else {
      taxa_table <- as.matrix(taxa_table)
    }
  }

  taxa_table <- .coerce_numeric_taxa_matrix(taxa_table)
  list(taxa_table = taxa_table, taxonomy_table = taxonomy_table)
}

.coerce_numeric_taxa_matrix <- function(x) {
  original_na <- is.na(x)
  if (!is.numeric(x)) {
    suppressWarnings(storage.mode(x) <- "double")
    if (any(is.na(x) & !original_na)) {
      stop("'taxa_table' must contain only numeric abundance values.")
    }
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("'taxa_table' must contain finite, non-missing abundance values.")
  }
  x
}

.align_taxa_table_to_dist <- function(taxa_table, dm) {
  dm_size <- attr(dm, "Size")
  dm_labels <- attr(dm, "Labels")
  sample_names_are_meaningful <- .has_meaningful_rownames(taxa_table)

  if (!is.null(dm_labels) && sample_names_are_meaningful) {
    if (anyDuplicated(rownames(taxa_table))) {
      stop("Sample names in 'taxa_table' must be unique.")
    }
    missing_samples <- setdiff(dm_labels, rownames(taxa_table))
    if (length(missing_samples) > 0) {
      stop("'taxa_table' is missing samples found in 'dm': ",
           paste(missing_samples, collapse = ", "), ".")
    }
    taxa_table <- taxa_table[dm_labels, , drop = FALSE]
  } else {
    if (nrow(taxa_table) != dm_size) {
      stop("'taxa_table' must contain the same number of samples as 'dm'.")
    }
    if (!is.null(dm_labels)) {
      rownames(taxa_table) <- dm_labels
    }
  }

  taxa_table
}

.name_taxa_columns <- function(taxa_table) {
  taxa_ids <- colnames(taxa_table)
  if (is.null(taxa_ids)) {
    taxa_ids <- paste0("Taxon_", seq_len(ncol(taxa_table)))
  }
  blank_taxa <- is.na(taxa_ids) | !nzchar(taxa_ids)
  taxa_ids[blank_taxa] <- paste0("Taxon_", which(blank_taxa))
  if (anyDuplicated(taxa_ids)) {
    stop("Taxon names must be unique.")
  }
  colnames(taxa_table) <- taxa_ids
  taxa_table
}

.select_taxonomy_output <- function(taxonomy_table, taxa_ids, taxonomy.ranks) {
  if (is.null(taxonomy.ranks) || length(taxonomy.ranks) == 0) {
    return(data.frame(row.names = seq_along(taxa_ids)))
  }
  if (!is.character(taxonomy.ranks) || anyNA(taxonomy.ranks)) {
    stop("'taxonomy.ranks' must be NULL, 'auto', 'all', or a character vector of rank names.")
  }
  if (is.null(taxonomy_table)) {
    if (identical(taxonomy.ranks, "auto")) {
      return(data.frame(row.names = seq_along(taxa_ids)))
    }
    stop("'taxonomy_table' is required when explicit 'taxonomy.ranks' are requested.")
  }

  taxonomy_rownames <- rownames(taxonomy_table)
  taxonomy_table <- data.frame(taxonomy_table, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(names(taxonomy_table))) {
    names(taxonomy_table) <- paste0("Rank_", seq_len(ncol(taxonomy_table)))
  }
  blank_ranks <- is.na(names(taxonomy_table)) | !nzchar(names(taxonomy_table))
  names(taxonomy_table)[blank_ranks] <- paste0("Rank_", which(blank_ranks))

  if (!is.null(taxonomy_rownames) &&
      !identical(taxonomy_rownames, as.character(seq_len(nrow(taxonomy_table))))) {
    rownames(taxonomy_table) <- taxonomy_rownames
  } else if (nrow(taxonomy_table) == length(taxa_ids)) {
    rownames(taxonomy_table) <- taxa_ids
  } else {
    stop("'taxonomy_table' must have row names matching taxa or the same number of rows as taxa.")
  }
  if (anyDuplicated(rownames(taxonomy_table))) {
    stop("Taxon names in 'taxonomy_table' must be unique.")
  }
  missing_taxa <- setdiff(taxa_ids, rownames(taxonomy_table))
  if (length(missing_taxa) > 0) {
    stop("'taxonomy_table' is missing taxa found in 'taxa_table': ",
         paste(missing_taxa, collapse = ", "), ".")
  }

  taxonomy_table <- taxonomy_table[match(taxa_ids, rownames(taxonomy_table)), , drop = FALSE]
  rownames(taxonomy_table) <- NULL

  selected_ranks <- .resolve_taxonomy_ranks(taxonomy_table, taxonomy.ranks)
  if (length(selected_ranks) == 0) {
    return(data.frame(row.names = seq_along(taxa_ids)))
  }

  taxonomy_table <- taxonomy_table[, selected_ranks, drop = FALSE]
  for (rank_name in names(taxonomy_table)) {
    taxonomy_table[[rank_name]] <- as.character(taxonomy_table[[rank_name]])
  }
  taxonomy_table
}

.resolve_taxonomy_ranks <- function(taxonomy_table, taxonomy.ranks) {
  if (length(taxonomy.ranks) == 1 && identical(taxonomy.ranks, "auto")) {
    rank_has_values <- vapply(
      taxonomy_table,
      function(x) any(!is.na(x) & nzchar(trimws(as.character(x)))),
      logical(1)
    )
    return(utils::tail(names(taxonomy_table)[rank_has_values], 2))
  }
  if (length(taxonomy.ranks) == 1 && identical(taxonomy.ranks, "all")) {
    return(names(taxonomy_table))
  }

  missing_ranks <- setdiff(taxonomy.ranks, names(taxonomy_table))
  if (length(missing_ranks) > 0) {
    stop("'taxonomy_table' is missing requested ranks: ",
         paste(missing_ranks, collapse = ", "), ".")
  }
  taxonomy.ranks
}

.summarize_taxa_matrix <- function(taxa_matrix) {
  data.frame(
    n.samples = rep(nrow(taxa_matrix), ncol(taxa_matrix)),
    taxon.mean = colMeans(taxa_matrix),
    taxon.prevalence = colMeans(taxa_matrix > 0),
    taxon.total.abundance = colSums(taxa_matrix),
    taxon.variance = vapply(
      seq_len(ncol(taxa_matrix)),
      function(i) stats::var(taxa_matrix[, i]),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

.validate_taxa_importance_formula <- function(formula) {
  if (is.null(formula)) {
    return(invisible(NULL))
  }
  # WdS.test() and a.dist() take right-hand-side formulas because the distance
  # matrix is supplied separately as dm, not as a left-hand-side response.
  if (!inherits(formula, "formula") || length(formula) != 2) {
    stop("'formula' must be NULL or a right-hand side formula such as ~ Subject_ID.")
  }
  invisible(NULL)
}

.prepare_taxa_importance_formula_data <- function(formula, formula_data, dm) {
  if (is.null(formula)) {
    # Taxon-only scans build a tiny formula_data object inside the loop, one
    # taxon at a time, so no sample metadata is needed here.
    return(NULL)
  }

  data <- .as_formula_data(formula_data)
  if (is.environment(data)) {
    # Environments support parent-frame lookup but do not have row names, so
    # alignment by sample ID is only possible for data-frame-like inputs.
    return(data)
  }

  .align_formula_data_to_dist(data, dm)
}

.align_formula_data_to_dist <- function(formula_data, dm) {
  dm_size <- attr(dm, "Size")
  dm_labels <- attr(dm, "Labels")

  if (!is.null(dm_labels) && .has_meaningful_rownames(formula_data)) {
    # Prefer label-based alignment when possible. This protects phyloseq-style
    # sample_data objects if their rows are in a different order than dm.
    if (anyDuplicated(rownames(formula_data))) {
      stop("Sample names in 'formula_data' must be unique.")
    }
    missing_samples <- setdiff(dm_labels, rownames(formula_data))
    if (length(missing_samples) > 0) {
      stop("'formula_data' is missing samples found in 'dm': ",
           paste(missing_samples, collapse = ", "), ".")
    }
    formula_data <- formula_data[dm_labels, , drop = FALSE]
  } else {
    # Without meaningful row names, fall back to positional matching but require
    # the same sample count as dm to avoid recycling or partial matching.
    if (nrow(formula_data) != dm_size) {
      stop("'formula_data' must contain the same number of samples as 'dm'.")
    }
    if (!is.null(dm_labels)) {
      rownames(formula_data) <- dm_labels
    }
  }

  formula_data
}

.build_taxa_importance_formula <- function(formula, formula_data,
                                           taxon_variable) {
  if (is.null(formula)) {
    # Backward-compatible default: one WdS.test() adjustment per taxon.
    return(stats::as.formula(paste("~", taxon_variable)))
  }
  if (taxon_variable %in% all.vars(formula)) {
    # Advanced users can place the taxon term exactly where they want it, for
    # example in an interaction or a no-intercept formula.
    return(formula)
  }

  # If the user supplies only sample-level adjustment terms, append the current
  # taxon term so the repeated scan still evaluates one taxon per WdS.test() fit.
  .add_factor_to_formula(formula, taxon_variable, formula_data)
}

.add_taxon_to_formula_data <- function(data, taxon_values, taxon_variable,
                                       sample_names) {
  if (is.null(data)) {
    data <- data.frame(taxon_values, check.names = FALSE)
    # Some vector/matrix inputs can carry the taxon ID into the data-frame
    # column name; force the reserved formula name so WdS.test() can always
    # resolve ~ .wdstar_taxon.
    names(data) <- taxon_variable
    rownames(data) <- sample_names
    return(data)
  }

  if (is.environment(data)) {
    # Use a child environment so each iteration can expose .wdstar_taxon without
    # modifying the user's environment or masking their original variables.
    taxon_env <- new.env(parent = data)
    assign(taxon_variable, taxon_values, envir = taxon_env)
    return(taxon_env)
  }

  # Data frames are copied in this local scope, so adding the placeholder column
  # here does not alter the user's original metadata object.
  data[[taxon_variable]] <- taxon_values
  data
}

.taxa_importance_row <- function(taxon, tested, skip.reason,
                                 adjustment.formula = NA_character_,
                                 adjustment.goodness.of.fit = NA_real_,
                                 factor.goodness.of.fit.unadjusted = NA_real_,
                                 factor.goodness.of.fit.after.adjustment = NA_real_,
                                 delta.factor.goodness.of.fit = NA_real_,
                                 full.goodness.of.fit = NA_real_,
                                 semi.partial.goodness.of.fit = NA_real_,
                                 partial.goodness.of.fit = NA_real_,
                                 statistic = NA_real_,
                                 p.value = NA_real_,
                                 omega.squared = NA_real_,
                                 omega.squared.unadjusted = NA_real_,
                                 delta.omega.squared = NA_real_,
                                 taxon_summary) {
  data.table::data.table(
    taxon = taxon,
    tested = tested,
    skip.reason = skip.reason,
    adjustment.formula = adjustment.formula,
    adjustment.goodness.of.fit = adjustment.goodness.of.fit,
    factor.goodness.of.fit.unadjusted = factor.goodness.of.fit.unadjusted,
    factor.goodness.of.fit.after.adjustment = factor.goodness.of.fit.after.adjustment,
    delta.factor.goodness.of.fit = delta.factor.goodness.of.fit,
    full.goodness.of.fit = full.goodness.of.fit,
    semi.partial.goodness.of.fit = semi.partial.goodness.of.fit,
    partial.goodness.of.fit = partial.goodness.of.fit,
    statistic = statistic,
    p.value = p.value,
    omega.squared = omega.squared,
    omega.squared.unadjusted = omega.squared.unadjusted,
    delta.omega.squared = delta.omega.squared,
    n.samples = taxon_summary$n.samples,
    taxon.mean = taxon_summary$taxon.mean,
    taxon.prevalence = taxon_summary$taxon.prevalence,
    taxon.total.abundance = taxon_summary$taxon.total.abundance,
    taxon.variance = taxon_summary$taxon.variance
  )
}

.run_taxa_importance_wdstest <- function(args, dots, verbose) {
  call_args <- c(args, dots)
  if (isTRUE(verbose)) {
    return(do.call(WdS.test, call_args))
  }
  suppressMessages(do.call(WdS.test, call_args))
}

.extract_goodness_value <- function(test_result, component) {
  components <- test_result$goodness.components
  if (is.null(components)) {
    return(NA_real_)
  }
  component_index <- match(component, components$component)
  if (is.na(component_index)) {
    return(NA_real_)
  }
  unname(components$pseudo.R2[[component_index]])
}

.validate_logical_scalar <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("'", arg, "' must be TRUE or FALSE.")
  }
}

.has_meaningful_rownames <- function(x) {
  row_names <- rownames(x)
  !is.null(row_names) && !identical(row_names, as.character(seq_len(nrow(x))))
}
