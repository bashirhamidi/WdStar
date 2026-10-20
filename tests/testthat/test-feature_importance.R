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

test_that("importance interfaces resolve metadata columns and external vectors", {
  metadata <- data.frame(
    treatment = factor(rep(c("A", "B"), 6)),
    block = factor(rep(1:3, each = 4)), age = 1:12,
    # Internal forwarding must not mistake these columns for resolved vectors.
    group = factor(rep("wrong", 12)), strata = rep("wrong", 12),
    row.names = paste0("sample", 1:12)
  )
  taxa <- cbind(ASV1 = sin(1:12), ASV2 = cos(1:12))
  rownames(taxa) <- rownames(metadata)
  dm <- dist(taxa)
  for (fun in list(WdS.feature.importance, WdS.taxa.importance)) {
    extra <- if (identical(fun, WdS.taxa.importance)) {
      list(taxa_table = taxa, taxa_are_rows = FALSE)
    } else list()
    run <- function(g, s, d) {
      set.seed(71)
      suppressMessages(do.call(fun, c(list(dm = dm, group = g, strata = s,
        adjustment = ~ age, data = d, nrep = 3, sort = FALSE), extra)))
    }
    expected <- run(metadata$treatment, metadata$block, metadata)
    expect_equal(run("treatment", "block", metadata)$importance, expected$importance)
    set.seed(71)
    bare <- suppressMessages(do.call(fun, c(list(dm = dm,
      group = quote(treatment), strata = quote(block), adjustment = ~ age,
      data = metadata, nrep = 3, sort = FALSE), extra)))
    expect_equal(bare$statistic, expected$statistic)
    expect_equal(bare$p.value, expected$p.value)
    if (identical(fun, WdS.taxa.importance)) {
      shuffled <- run("treatment", "block", metadata[12:1, ])
      expect_equal(shuffled$importance, expected$importance)
      expect_equal(shuffled$p.value, expected$p.value)
    }
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      from_phyloseq <- run("treatment", "block", phyloseq::sample_data(metadata))
      expect_equal(from_phyloseq$importance, expected$importance)
    }
    # NULL metadata uses the caller, including for adjustment variables.
    age <- metadata$age
    expect_equal(run(metadata$treatment, metadata$block, NULL)$importance,
                 expected$importance)
  }
})

test_that("single-feature models honor diagnostic options on adjusted inputs", {
  covariate <- seq_len(fi_n_rows)
  adjusted_dm <- suppressMessages(a.dist(fi_distance_matrix, ~ covariate))
  result <- suppressMessages(WdS.feature.importance(
    adjusted_dm, fi_factor_var, ~ covariate, nrep = 1,
    distance.diagnostics = FALSE
  ))
  expect_null(result$test.result[[1]]$distance.diagnostics)
})

test_that("importance interfaces reject old named arguments before fitting", {
  for (fun in list(WdS.feature.importance, WdS.taxa.importance)) {
    for (old_name in c("f", "formula", "formula_data")) {
      args <- list(dm = fi_distance_matrix, group = fi_factor_var)
      args[old_name] <- list(NULL)
      expect_error(do.call(fun, args), "Unsupported argument", fixed = TRUE)
    }
    expect_error(fun(dm = fi_distance_matrix, group = fi_factor_var,
                     goodness = "all"), "controlled internally")
  }
})

test_that("WdS.feature.importance ranks leave-one-feature-out results", {
  importance_data <- data.frame(
    adjustment1 = rep(c("low", "high"), length.out = fi_n_rows),
    adjustment2 = seq_len(fi_n_rows),
    adjustment3 = rep(c("yes", "no", "maybe", "yes"), length.out = fi_n_rows)
  )

  result <- suppressMessages(WdS.feature.importance(
    dm = fi_distance_matrix,
    group = fi_factor_var,
    adjustment = ~ adjustment1 + adjustment2 + adjustment3,
    data = importance_data,
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
    group = fi_factor_var,
    adjustment = ~ adjustment,
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
    WdS.feature.importance(dm = fi_distance_matrix, group = fi_factor_var, nrep = 9),
    "'adjustment' must be provided"
  )
  expect_error(
    WdS.feature.importance(
      dm = fi_distance_matrix,
      group = fi_factor_var,
      adjustment = ~ 1,
      nrep = 9
    ),
    "at least one adjustment feature"
  )
})

test_that("WdS.taxa.importance ranks taxa and includes automatic taxonomy ranks", {
  taxa_table <- data.frame(
    ASV1 = fi_example_df$X1,
    ASV2 = fi_example_df$X2,
    ASV3 = fi_example_df$X3,
    ASV4 = 1,
    row.names = attr(fi_distance_matrix, "Labels")
  )
  taxonomy_table <- data.frame(
    Kingdom = rep("Bacteria", 4),
    Phylum = c("Firmicutes", "Bacteroidota", "Proteobacteria", "Actinobacteriota"),
    Genus = c("GenusA", "GenusB", "GenusC", "GenusD"),
    Species = c("species_a", "species_b", "species_c", "species_d"),
    row.names = colnames(taxa_table)
  )

  result <- WdS.taxa.importance(
    dm = fi_distance_matrix,
    group = fi_factor_var,
    taxa_table = taxa_table,
    taxa_are_rows = FALSE,
    taxonomy_table = taxonomy_table,
    nrep = 9
  )

  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 4)
  expect_true(all(c("taxon", "Genus", "Species", "importance") %in% names(result)))
  expect_false("Kingdom" %in% names(result))
  expect_false("Phylum" %in% names(result))
  expect_false("test.result" %in% names(result))
  expect_equal(
    result$importance,
    result$delta.factor.goodness.of.fit
  )
  expect_true(all(diff(result$importance[!is.na(result$importance)]) <= 0))

  constant_taxon <- result[result$taxon == "ASV4", ]
  expect_false(constant_taxon$tested)
  expect_equal(constant_taxon$skip.reason, "zero variance")
  expect_true(is.na(constant_taxon$importance))
})

test_that("WdS.taxa.importance accepts selected taxonomy ranks and keeps tests", {
  taxa_table <- t(as.matrix(data.frame(
    ASV1 = fi_example_df$X1,
    ASV2 = fi_example_df$X2,
    row.names = attr(fi_distance_matrix, "Labels")
  )))
  taxonomy_table <- data.frame(
    Family = c("FamilyA", "FamilyB"),
    Genus = c("GenusA", "GenusB"),
    Species = c("species_a", "species_b"),
    row.names = rownames(taxa_table)
  )

  result <- WdS.taxa.importance(
    dm = fi_distance_matrix,
    group = fi_factor_var,
    taxa_table = taxa_table,
    taxa_are_rows = TRUE,
    taxonomy_table = taxonomy_table,
    taxonomy.ranks = c("Family", "Genus"),
    rank.by = "adjustment.goodness.of.fit",
    keep.tests = TRUE,
    sort = FALSE,
    nrep = 9
  )

  expect_s3_class(result, "data.table")
  expect_equal(result$taxon, c("ASV1", "ASV2"))
  expect_true(all(c("Family", "Genus", "test.result") %in% names(result)))
  expect_false("Species" %in% names(result))
  expect_true(all(result$rank.by == "adjustment.goodness.of.fit"))
  expect_equal(result$importance, result$adjustment.goodness.of.fit)
  expect_true(all(vapply(result$test.result, inherits, logical(1), what = "wdstest")))
})

test_that("WdS.taxa.importance adds formula terms to each taxon adjustment", {
  taxa_table <- data.frame(
    ASV1 = fi_example_df$X1,
    ASV2 = fi_example_df$X2,
    row.names = attr(fi_distance_matrix, "Labels")
  )
  sample_metadata <- data.frame(
    Subject_ID = factor(rep(seq_len(fi_n_rows / 2), each = 2)),
    row.names = rev(attr(fi_distance_matrix, "Labels"))
  )
  # Reverse both row names and values so the test proves data is aligned
  # by sample ID before the current taxon is added.
  sample_metadata$Subject_ID <- rev(sample_metadata$Subject_ID)

  result <- WdS.taxa.importance(
    dm = fi_distance_matrix,
    group = fi_factor_var,
    taxa_table = taxa_table,
    taxa_are_rows = FALSE,
    adjustment = ~ Subject_ID,
    data = sample_metadata,
    rank.by = "adjustment.goodness.of.fit",
    sort = FALSE,
    nrep = 9
  )

  manual_data <- data.frame(
    ASV = taxa_table$ASV1,
    Subject_ID = sample_metadata[rownames(taxa_table), "Subject_ID"],
    row.names = rownames(taxa_table)
  )
  manual_result <- suppressMessages(WdS.test(
    dm = fi_distance_matrix,
    group = fi_factor_var,
    # WdS.taxa.importance(adjustment = ~ Subject_ID) should be equivalent to this
    # manual model after it adds the current taxon to the adjustment formula.
    adjustment = ~ Subject_ID + ASV,
    data = manual_data,
    goodness = c("adjustment", "semi.partial"),
    nrep = 9
  ))

  expect_equal(
    result[result[["taxon"]] == "ASV1", ][["adjustment.goodness.of.fit"]],
    .extract_goodness_value(manual_result, "adjustment")
  )
  expect_equal(
    result[result[["taxon"]] == "ASV1", ][["semi.partial.goodness.of.fit"]],
    .extract_goodness_value(manual_result, "semi.partial")
  )
  expect_match(result$adjustment.formula[[1]], "\\.wdstar_taxon")
  expect_match(result$adjustment.formula[[1]], "Subject_ID")
})

test_that("WdS.taxa.importance accepts explicit taxon placeholders in formula", {
  taxa_table <- data.frame(
    ASV1 = fi_example_df$X1,
    row.names = attr(fi_distance_matrix, "Labels")
  )
  sample_metadata <- data.frame(
    Subject_ID = factor(rep(seq_len(fi_n_rows / 2), each = 2)),
    row.names = attr(fi_distance_matrix, "Labels")
  )

  result <- WdS.taxa.importance(
    dm = fi_distance_matrix,
    group = fi_factor_var,
    taxa_table = taxa_table,
    taxa_are_rows = FALSE,
    # Including .wdstar_taxon explicitly lets users control the formula shape
    # instead of having the taxon term appended automatically.
    adjustment = ~ .wdstar_taxon + Subject_ID,
    data = sample_metadata,
    sort = FALSE,
    nrep = 9
  )

  manual_data <- data.frame(
    ASV = taxa_table$ASV1,
    Subject_ID = sample_metadata$Subject_ID,
    row.names = rownames(taxa_table)
  )
  manual_result <- suppressMessages(WdS.test(
    dm = fi_distance_matrix,
    group = fi_factor_var,
    adjustment = ~ ASV + Subject_ID,
    data = manual_data,
    goodness = c("adjustment", "semi.partial"),
    nrep = 9
  ))

  expect_equal(
    result$adjustment.goodness.of.fit,
    .extract_goodness_value(manual_result, "adjustment")
  )
  expect_equal(
    result$semi.partial.goodness.of.fit,
    .extract_goodness_value(manual_result, "semi.partial")
  )
})

test_that("WdS.taxa.importance validates taxa and taxonomy inputs", {
  expect_error(
    WdS.taxa.importance(dm = fi_distance_matrix, group = fi_factor_var, nrep = 9),
    "Provide either 'physeq' or 'taxa_table'"
  )

  bad_taxa_table <- data.frame(
    ASV1 = seq_len(fi_n_rows - 1)
  )
  expect_error(
    WdS.taxa.importance(
      dm = fi_distance_matrix,
      group = fi_factor_var,
      taxa_table = bad_taxa_table,
      taxa_are_rows = FALSE,
      nrep = 9
    ),
    "same number of samples"
  )

  taxa_table <- data.frame(
    ASV1 = fi_example_df$X1,
    ASV2 = fi_example_df$X2,
    row.names = attr(fi_distance_matrix, "Labels")
  )
  taxonomy_table <- data.frame(
    Genus = "GenusA",
    row.names = "ASV1"
  )
  expect_error(
    WdS.taxa.importance(
      dm = fi_distance_matrix,
      group = fi_factor_var,
      taxa_table = taxa_table,
      taxa_are_rows = FALSE,
      taxonomy_table = taxonomy_table,
      taxonomy.ranks = "Genus",
      nrep = 9
    ),
    "missing taxa"
  )

  expect_error(
    WdS.taxa.importance(
      dm = fi_distance_matrix,
      group = fi_factor_var,
      taxa_table = taxa_table,
      taxa_are_rows = FALSE,
      adjustment = Condition ~ Subject_ID,
      data = data.frame(
        Subject_ID = seq_len(fi_n_rows),
        Condition = fi_factor_var
      ),
      nrep = 9
    ),
    "right-hand-side formula"
  )

  expect_error(
    WdS.taxa.importance(
      dm = fi_distance_matrix,
      group = fi_factor_var,
      taxa_table = taxa_table,
      taxa_are_rows = FALSE,
      adjustment = ~ Subject_ID,
      data = data.frame(Subject_ID = seq_len(fi_n_rows - 1)),
      nrep = 9
    ),
    "same number of samples"
  )
})
