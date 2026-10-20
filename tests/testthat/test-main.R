test_that("WdS accepts prepared factors and agrees with WdS.test", {
  set.seed(92)
  dm <- dist(matrix(rnorm(48), nrow = 12))
  metadata <- data.frame(treatment = factor(rep(c("A", "B"), each = 6)))
  group_vector <- metadata$treatment
  expected <- WdS(dm, group_vector)
  expect_equal(WdS(dm, group = group_vector), expected)
  expect_warning(legacy <- WdS(dm, f = group_vector), "'f' is deprecated in WdS()",
                 fixed = TRUE)
  expect_equal(legacy, expected)
  # The permutation engine still calls the statistic positionally.
  result <- WdS.test(dm, treatment, data = metadata, nrep = 3, goodness = "none")
  expect_equal(unname(result$statistic), unname(expected))
})

test_that("WdS rejects ambiguous and invalid grouping arguments", {
  dm <- dist(matrix(seq_len(24), nrow = 12))
  group <- factor(rep(1:2, each = 6))
  expect_error(WdS(dm), "'group' must be provided", fixed = TRUE)
  expect_error(WdS(dm, group, f = group), "Use only one")
  expect_error(WdS(dm, f = group, f = group), "more than once")
  expect_error(WdS(dm, group, typo = TRUE), "Only")
  expect_error(WdS(dm, group, data = data.frame(group)), "Only")
  expect_error(WdS(dm, as.character(group)), "must be a factor")
  expect_error(WdS(dm, group[-1]), "must have length 12")
  expect_error(WdS(dm, rep("A", 12)), "at least two levels")
})
