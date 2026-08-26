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
#' data(mtcars)
#' dm <- dist(mtcars[1:3], method = "euclidean")
#' f <- factor(mtcars$gear)
#' WdS.feature.importance(
#'   dm = dm,
#'   f = f,
#'   formula = ~ wt + as.factor(am),
#'   formula_data = mtcars,
#'   nrep = 9
#' )
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
