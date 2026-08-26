set.seed(44)
n_rows <- 100

example_df <- data.frame(
  X1 = rnorm(n_rows, mean = 0, sd = 1),
  X2 = rnorm(n_rows, mean = 0, sd = 1),
  X3 = rnorm(n_rows, mean = 0, sd = 1),
  X4 = rnorm(n_rows, mean = 0, sd = 1)
)

# Generate factors
factor_var <- factor(rep(c("A", "B"), each = n_rows / 2))

# Make distance matrix (type: dist)
distance_matrix <- dist(example_df)

# Make square matrix form (type: matrix, array)
distance_matrix_sq <- as.matrix(distance_matrix)

ss2 <- dist.ss2(distance_matrix_sq, factor_var)
cat("Sum of Squares Matrix:\n")
print(ss2)

test_that("'distance' objects are successfully identified", {
  expect_true(is.dist(distance_matrix))
  expect_false(is.dist(distance_matrix_sq))
})

test_that("correct sigma squared value is calculated", {
  expect_equal(dist.sigma2(distance_matrix), 4.2206406)
})

test_that("dist.ss2 correctly calcualtes sum of squares matrix", {
  expect_equal(dist.ss2(distance_matrix_sq, factor_var)[1], 3323.9046)
})

test_that("dist.group.sigma2 correctly calcualtes grou-wise sigma squared stats", {
  expect_equal(dist.group.sigma2(distance_matrix, factor_var)[1], setNames(4.1662956, "A"))
  expect_equal(dist.group.sigma2(distance_matrix, factor_var)[2], setNames(4.2468826, "B"))
})

test_that("dist.cohen.d correctly calculates Cohen's d for the distance matrix", {
  expect_equal(dist.cohen.d(distance_matrix, factor_var), setNames(0.23071167, "A"),
               tolerance = 0.001)
})

test_that("a.dist uses formula variables from the parent frame by default", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment))

  expect_s3_class(adjusted_dm, "dist")
  expect_equal(attr(adjusted_dm, "Size"), n_rows)
})

test_that("WdS.test uses formula variables from the parent frame by default", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment
  ))

  expect_s3_class(result, "htest")
  expect_match(result$method, "Adjusted")
})

test_that("WdS.test accepts phyloseq sample_data as formula_data", {
  skip_if_not_installed("phyloseq")

  sample_df <- data.frame(
    adjustment = rep(c("low", "high"), length.out = n_rows),
    row.names = attr(distance_matrix, "Labels")
  )
  sample_df <- phyloseq::sample_data(sample_df)

  result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    formula_data = sample_df
  ))

  expect_s3_class(result, "htest")
  expect_match(result$method, "Adjusted")
})

test_that("WdS.test reports auto goodness of fit for adjusted and unadjusted tests", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment))
  factor_data <- data.frame(.wdstar_f = factor_var)
  factor_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ .wdstar_f, formula_data = factor_data))
  full_data <- data.frame(adjustment = adjustment, .wdstar_f = factor_var)
  full_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment + .wdstar_f, formula_data = full_data))
  expected_factor_r_squared <- dist.goodness.of.fit(distance_matrix, factor_dm)
  expected_adjustment_r_squared <- dist.goodness.of.fit(distance_matrix, adjusted_dm)
  expected_full_r_squared <- dist.goodness.of.fit(distance_matrix, full_dm)
  expected_semi_partial <- unname(expected_full_r_squared) - unname(expected_adjustment_r_squared)
  expected_partial <- expected_semi_partial / (1 - unname(expected_adjustment_r_squared))

  unadjusted_result <- WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9
  )
  adjusted_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment
  ))

  expect_named(unadjusted_result$estimate, "effect size estimator of variance, omega-squared (\u03C9\u00B2)")
  expect_named(adjusted_result$estimate, "effect size estimator of variance, omega-squared (\u03C9\u00B2)")
  expect_named(unadjusted_result$goodness.of.fit, "distance-based pseudo-R\u00B2 for factor only")
  expect_named(adjusted_result$goodness.of.fit, "semi-partial distance-based pseudo-R\u00B2 for factor given adjustment")
  expect_equal(unname(unadjusted_result$goodness.of.fit), unname(expected_factor_r_squared))
  expect_equal(unname(adjusted_result$goodness.of.fit), expected_semi_partial)
  expect_setequal(
    adjusted_result$goodness.components$component,
    c("adjustment", "full", "semi.partial", "partial")
  )
  expect_equal(
    adjusted_result$goodness.components$pseudo.R2[
      adjusted_result$goodness.components$component == "adjustment"
    ],
    unname(expected_adjustment_r_squared)
  )
  expect_equal(
    adjusted_result$goodness.components$pseudo.R2[
      adjusted_result$goodness.components$component == "full"
    ],
    unname(expected_full_r_squared)
  )
  expect_equal(
    adjusted_result$goodness.components$pseudo.R2[
      adjusted_result$goodness.components$component == "partial"
    ],
    expected_partial
  )
  expect_equal(
    adjusted_result$goodness.components$component[adjusted_result$goodness.components$displayed],
    "semi.partial"
  )
  expect_equal(
    adjusted_result$goodness.components$formula[
      adjusted_result$goodness.components$component == "full"
    ],
    "~ adjustment + f"
  )
  expect_equal(adjusted_result$distance.diagnostics$role, c("adjustment", "full"))
  expect_false("eigenvalues" %in% names(adjusted_result$distance.diagnostics))
})

test_that("WdS.test can report selected goodness-of-fit components", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  adjustment_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    goodness = "adjustment"
  ))
  all_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    goodness = "all"
  ))
  none_result <- WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    goodness = "none"
  )

  expect_named(adjustment_result$goodness.of.fit, "distance-based pseudo-R\u00B2 for adjustment variables")
  expect_setequal(
    all_result$goodness.components$component,
    c("adjustment", "factor.only", "full", "semi.partial", "partial")
  )
  expect_equal(
    all_result$goodness.components$component[all_result$goodness.components$displayed],
    c("adjustment", "factor.only", "full", "semi.partial", "partial")
  )
  expect_null(none_result$goodness.of.fit)
  expect_null(none_result$goodness.components)
  expect_null(none_result$distance.diagnostics)
  expect_named(
    WdS.test(
      dm = distance_matrix,
      f = factor_var,
      nrep = 9,
      goodness = stats::setNames("auto", "user_supplied_name")
    )$goodness.of.fit,
    "distance-based pseudo-R\u00B2 for factor only"
  )
  expect_error(
    WdS.test(dm = distance_matrix, f = factor_var, nrep = 9, goodness = "adjustment"),
    "requires an adjustment 'formula'"
  )
  expect_error(
    WdS.test(dm = distance_matrix, f = factor_var, nrep = 9, goodness = NA_character_),
    "cannot contain missing values"
  )
})

test_that("WdS.test collects distance diagnostics separately from goodness components", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  pre_adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment))

  adjustment_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    goodness = "adjustment"
  ))
  all_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    goodness = "all",
    keep.eigenvalues = TRUE
  ))
  no_diagnostics_result <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    distance.diagnostics = FALSE
  ))
  pre_adjusted_result <- WdS.test(
    dm = pre_adjusted_dm,
    f = factor_var,
    nrep = 9,
    goodness = "none"
  )

  expect_equal(adjustment_result$distance.diagnostics$role, "adjustment")
  expect_equal(
    all_result$distance.diagnostics$role,
    c("adjustment", "factor.only", "full")
  )
  expect_true("eigenvalues" %in% names(all_result$distance.diagnostics))
  expect_length(all_result$distance.diagnostics$eigenvalues[[1]], n_rows)
  expect_null(no_diagnostics_result$distance.diagnostics)
  expect_equal(pre_adjusted_result$distance.diagnostics$role, "input")
  expect_equal(pre_adjusted_result$distance.diagnostics$formula, "~adjustment")
})

test_that("WdS.test goodness formula construction preserves no-intercept formulas", {
  formula_data <- data.frame(adjustment = seq_len(n_rows))
  full_formula <- WdStar:::.add_factor_to_formula(~ adjustment - 1, ".wdstar_f", formula_data)

  expect_equal(deparse(full_formula), "~adjustment + .wdstar_f - 1")
  expect_equal(
    WdStar:::.display_formula_with_factor(~ adjustment - 1, formula_data),
    "~ adjustment + f - 1"
  )
})

test_that("WdS.test suppresses internal goodness residualization messages", {
  expect_silent(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9
  ))
})

test_that("estimate labels fall back safely for non-UTF-8 HTML output", {
  x <- WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9
  )

  expect_match(names(x$estimate), "\u03C9\u00B2", fixed = TRUE)
  expect_equal(
    WdStar:::.format_estimate_name(
      names(x$estimate), utf8 = FALSE, html = TRUE, latex = FALSE
    ),
    "effect size estimator of variance, omega-squared (omega^2)"
  )
  expect_equal(
    WdStar:::.format_estimate_name(
      names(x$estimate), utf8 = TRUE, html = TRUE, latex = FALSE
    ),
    names(x$estimate)
  )
  expect_equal(
    WdStar:::.format_estimate_name(
      names(x$estimate), utf8 = TRUE, html = FALSE, latex = TRUE
    ),
    "effect size estimator of variance, omega-squared (omega^2)"
  )
})

test_that("dist.goodness.of.fit computes distance-based pseudo-R squared", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment))
  expected_r_squared <- 1 - (dist.sigma2(adjusted_dm) / dist.sigma2(distance_matrix))

  result <- dist.goodness.of.fit(distance_matrix, dm_residual = adjusted_dm)
  compatibility_result <- dist.goodness.of.fit(distance_matrix, adjusted_dm = adjusted_dm)

  expect_named(result, "distance-based pseudo-R\u00B2")
  expect_equal(unname(result), expected_r_squared)
  expect_equal(unname(compatibility_result), expected_r_squared)
})

test_that("distance-based pseudo-R squared has expected boundary cases", {
  expect_equal(unname(dist.goodness.of.fit(distance_matrix, distance_matrix)), 0)
  expect_equal(
    unname(dist.goodness.of.fit(distance_matrix, distance_matrix / 2)),
    0.75
  )

  zero_dm <- stats::as.dist(matrix(0, nrow = 4, ncol = 4))
  expect_true(is.na(unname(dist.goodness.of.fit(zero_dm, zero_dm))))
})

test_that("dist.goodness.of.fit validates distance matrix inputs", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment))

  expect_error(
    dist.goodness.of.fit(as.matrix(distance_matrix), adjusted_dm),
    "must both be distance matrices"
  )
  expect_error(
    dist.goodness.of.fit(distance_matrix, as.dist(matrix(0, nrow = 3, ncol = 3))),
    "same number of observations"
  )
  expect_error(
    dist.goodness.of.fit(distance_matrix),
    "'dm_residual' must be provided"
  )
})

test_that("a.dist stores compact distance diagnostics by default", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment, tol = 0.1))
  diagnostics <- attr(adjusted_dm, "distance.diagnostics")

  expect_s3_class(diagnostics, "data.frame")
  expect_equal(nrow(diagnostics), 1)
  expect_equal(diagnostics$tol, 0.1)
  expect_equal(diagnostics$n.eigenvalues, n_rows)
  expect_gt(diagnostics$n.zeroed.by.tol, 0)
  expect_equal(
    diagnostics$n.negative.removed,
    diagnostics$n.negative.raw - diagnostics$n.negative.zeroed.by.tol
  )
  expect_true(all(c(
    "sum.positive.raw",
    "sum.negative.raw",
    "sum.negative.magnitude.raw",
    "prop.negative.raw"
  ) %in% names(diagnostics)))
  expect_false("eigenvalues" %in% names(diagnostics))
})

test_that("a.dist can keep raw eigenvalues when requested", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(
    distance_matrix,
    formula = ~ adjustment,
    keep.eigenvalues = TRUE
  ))
  diagnostics <- attr(adjusted_dm, "distance.diagnostics")

  expect_true(diagnostics$keep.eigenvalues)
  expect_true("eigenvalues" %in% names(diagnostics))
  expect_type(diagnostics$eigenvalues[[1]], "double")
  expect_length(diagnostics$eigenvalues[[1]], n_rows)
})

test_that("a.dist can skip distance diagnostics", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(
    distance_matrix,
    formula = ~ adjustment,
    distance.diagnostics = FALSE
  ))

  expect_null(attr(adjusted_dm, "distance.diagnostics", exact = TRUE))
  expect_error(
    a.dist(
      distance_matrix,
      formula = ~ adjustment,
      distance.diagnostics = FALSE,
      keep.eigenvalues = TRUE
    ),
    "requires 'distance.diagnostics = TRUE'"
  )
})

test_that("a.dist reports and zeroes eigenvalues below tolerance", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)

  expect_message(
    adjusted_dm <- a.dist(distance_matrix, formula = ~ adjustment, tol = 0.1),
    "eigenvalues were smaller than 'tol'"
  )

  expect_s3_class(adjusted_dm, "dist")
  expect_false(any(is.nan(adjusted_dm)))
})

test_that("a.dist reports and zeroes negative eigenvalues", {
  non_euclidean_dm <- as.dist(matrix(
    c(
      0.000, 2.114, 1.444, 1.293,
      2.114, 0.000, 1.133, 1.278,
      1.444, 1.133, 0.000, 2.059,
      1.293, 1.278, 2.059, 0.000
    ),
    nrow = 4
  ))
  attr(non_euclidean_dm, "method") <- "non-euclidean"
  adjustment <- c(1, 2, 1, 2)

  expect_message(
    adjusted_dm <- a.dist(non_euclidean_dm, formula = ~ adjustment, tol = 0),
    "eigenvalues were negative and were set to zero"
  )

  expect_s3_class(adjusted_dm, "dist")
  expect_false(any(is.nan(adjusted_dm)))
})
