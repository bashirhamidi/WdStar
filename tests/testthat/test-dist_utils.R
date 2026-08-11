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

test_that("WdS.test reports goodness of fit only for adjusted tests", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment))
  expected_r_squared <- dist.goodness.of.fit(distance_matrix, adjusted_dm)

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

  expect_null(unadjusted_result$goodness.of.fit)
  expect_named(unadjusted_result$estimate, "effect size estimator of variance, omega-squared (\u03C9\u00B2)")
  expect_named(adjusted_result$estimate, "effect size estimator of variance, omega-squared (\u03C9\u00B2)")
  expect_named(adjusted_result$goodness.of.fit, "goodness-of-fit coefficient of determination (R\u00B2)")
  expect_equal(unname(adjusted_result$goodness.of.fit), unname(expected_r_squared))
})

test_that("dist.goodness.of.fit computes selected goodness-of-fit metrics", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_dm <- suppressMessages(a.dist(distance_matrix, formula = ~ adjustment))
  expected_r_squared <- 1 - (dist.sigma2(adjusted_dm) / dist.sigma2(distance_matrix))
  residual_distances <- as.numeric(adjusted_dm)
  n <- length(residual_distances)
  rss <- sum(residual_distances^2)
  expected_aic <- 2 * 3 - 2 * (-0.5 * n * (log(2 * pi * (rss / n)) + 1))

  result <- dist.goodness.of.fit(distance_matrix, adjusted_dm)
  aic_result <- dist.goodness.of.fit(distance_matrix, adjusted_dm, metric = "aic", k = 3)
  both_result <- dist.goodness.of.fit(distance_matrix, adjusted_dm, metric = "both", k = 3)

  expect_named(result, "goodness-of-fit coefficient of determination (R\u00B2)")
  expect_equal(unname(result), expected_r_squared)
  expect_named(aic_result, "goodness-of-fit Akaike information criterion (AIC)")
  expect_equal(unname(aic_result), expected_aic)
  expect_equal(unname(both_result[1]), expected_r_squared)
  expect_equal(unname(both_result[2]), expected_aic)
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
    dist.goodness.of.fit(distance_matrix, adjusted_dm, metric = "aic", k = 0),
    "'k' must be a single positive numeric value"
  )
})

test_that("WdS.test supports selecting AIC and both goodness-of-fit metrics", {
  adjustment <- rep(c("low", "high"), length.out = n_rows)
  adjusted_result_aic <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    fit_metric = "aic"
  ))
  adjusted_result_both <- suppressMessages(WdS.test(
    dm = distance_matrix,
    f = factor_var,
    nrep = 9,
    formula = ~ adjustment,
    fit_metric = "both"
  ))

  expect_named(adjusted_result_aic$goodness.of.fit, "goodness-of-fit Akaike information criterion (AIC)")
  expect_named(
    adjusted_result_both$goodness.of.fit,
    c(
      "goodness-of-fit coefficient of determination (R\u00B2)",
      "goodness-of-fit Akaike information criterion (AIC)"
    )
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
