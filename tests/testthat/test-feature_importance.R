set.seed(44)
fi_n_rows <- 100
fi_example_df <- data.frame(
  X1 = rnorm(fi_n_rows, mean = 0, sd = 1),
  X2 = rnorm(fi_n_rows, mean = 0, sd = 1),
  X3 = rnorm(fi_n_rows, mean = 0, sd = 1),
  X4 = rnorm(fi_n_rows, mean = 0, sd = 1)
)
fi_factor_var <- factor(rep(c("A", "B"), each = fi_n_rows / 2))
fi_distance_matrix <- dist(fi_example_df)

test_that("WdS.feature.importance ranks leave-one-feature-out results", {
  importance_data <- data.frame(
    adjustment1 = rep(c("low", "high"), length.out = fi_n_rows),
    adjustment2 = seq_len(fi_n_rows),
    adjustment3 = rep(c("yes", "no", "maybe", "yes"), length.out = fi_n_rows)
  )

  result <- suppressMessages(WdS.feature.importance(
    dm = fi_distance_matrix,
    f = fi_factor_var,
    formula = ~ adjustment1 + adjustment2 + adjustment3,
    formula_data = importance_data,
    nrep = 9
  ))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_setequal(result$feature, c("adjustment1", "adjustment2", "adjustment3"))
  expect_true(all(diff(result$importance) <= 0))
  expect_equal(
    result$importance,
    result$full.goodness.of.fit - result$leave.one.out.goodness.of.fit
  )
  expect_true(all(vapply(result$test.result, inherits, logical(1), what = "wdstest")))
})

test_that("WdS.feature.importance handles formulas with one feature", {
  adjustment <- rep(c("low", "high"), length.out = fi_n_rows)

  result <- suppressMessages(WdS.feature.importance(
    dm = fi_distance_matrix,
    f = fi_factor_var,
    formula = ~ adjustment,
    nrep = 9
  ))

  expect_equal(nrow(result), 1)
  expect_equal(result$feature, "adjustment")
  expect_equal(result$leave.one.out.goodness.of.fit, 0)
  expect_equal(result$importance, result$full.goodness.of.fit)
  expect_true(is.na(result$formula))
  expect_null(result$test.result[[1]]$goodness.of.fit)
})

test_that("WdS.feature.importance validates formula input", {
  expect_error(
    WdS.feature.importance(dm = fi_distance_matrix, f = fi_factor_var, nrep = 9),
    "'formula' must be provided"
  )
  expect_error(
    WdS.feature.importance(
      dm = fi_distance_matrix,
      f = fi_factor_var,
      formula = ~ 1,
      nrep = 9
    ),
    "at least one adjustment feature"
  )
})
