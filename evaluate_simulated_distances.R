library(WdStar)
library(phyloseq)

simulate_microbiome <- function(
  seed = 1103,
  n = 120L,
  taxa = 100L,
  discordant_per_z_level = 8L,
  sequencing_depth = 5000L
) {
  set.seed(seed)
  z01 <- rep(c(0, 1), each = n / 2L)
  f01 <- z01
  for (z_value in 0:1) {
    candidates <- which(z01 == z_value)
    flip <- sample(candidates, discordant_per_z_level)
    f01[flip] <- 1L - f01[flip]
  }
  z_centered <- z01 - mean(z01)
  f_centered <- f01 - mean(f01)
  baseline <- rnorm(taxa, mean = -2.5, sd = 0.8)
  beta_z <- beta_f <- numeric(taxa)
  beta_z[1:15] <- 1
  beta_z[16:30] <- -1
  beta_f[1:15] <- -2
  beta_f[16:30] <- 2
  beta_f[31:45] <- 1
  beta_f[46:60] <- -1
  counts <- matrix(0L, nrow = taxa, ncol = n)
  for (i in seq_len(n)) {
    log_abundance <- baseline + z_centered[i] * beta_z +
      f_centered[i] * beta_f + rnorm(taxa, sd = 0.4)
    probabilities <- exp(log_abundance - max(log_abundance))
    probabilities <- probabilities / sum(probabilities)
    counts[, i] <- as.vector(rmultinom(1L, sequencing_depth, probabilities))
  }
  rownames(counts) <- paste0("Taxon", seq_len(taxa))
  colnames(counts) <- paste0("Sample", seq_len(n))
  metadata <- data.frame(
    Z = factor(z01, labels = c("Z0", "Z1")),
    F = factor(f01, labels = c("F0", "F1")),
    row.names = colnames(counts)
  )
  phyloseq(otu_table(counts, taxa_are_rows = TRUE), sample_data(metadata))
}

negative_fraction <- function(dm, tol = 1e-8) {
  d2 <- as.matrix(dm)^2
  n <- nrow(d2)
  J <- diag(n) - matrix(1 / n, n, n)
  G <- -0.5 * J %*% d2 %*% J
  values <- eigen((G + t(G)) / 2, symmetric = TRUE, only.values = TRUE)$values
  threshold <- tol * max(1, max(abs(values)))
  values <- values[abs(values) > threshold]
  sum(abs(values[values < 0])) / sum(abs(values))
}

evaluate_distance <- function(label, dm, metadata) {
  z_resid <- a.dist(dm, ~ Z, metadata)
  full_resid <- a.dist(dm, ~ Z + F, metadata)
  two_resid <- a.dist(z_resid, ~ F, metadata)
  ad <- vegan::adonis2(dm ~ Z + F, metadata, by = "terms", permutations = 0)
  r2_z <- as.numeric(dist.goodness.of.fit(dm, z_resid))
  one <- as.numeric(dist.goodness.of.fit(dm, full_resid))
  two <- as.numeric(dist.goodness.of.fit(dm, two_resid))
  vegan_full <- 1 - ad["Residual", "R2"]
  data.frame(
    distance = label,
    negative_absolute_fraction = negative_fraction(dm),
    vegan_Z = ad["Z", "R2"],
    vegan_full = vegan_full,
    vegan_F_after_Z = ad["F", "R2"],
    WdStar_Z = r2_z,
    WdStar_one_pass_full = one,
    WdStar_two_pass_full = two,
    one_minus_vegan = one - vegan_full,
    one_minus_two = one - two
  )
}

ps <- simulate_microbiome()
md <- data.frame(sample_data(ps))
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
proportions <- otu / rowSums(otu)

d_bray <- phyloseq::distance(ps, "bray")
d_jaccard_quantitative <- phyloseq::distance(ps, "jaccard")
d_jaccard_binary <- vegan::vegdist(otu, method = "jaccard", binary = TRUE)
d_js_divergence <- phyloseq::distance(ps, "jsd")
d_sqrt_js <- sqrt(d_js_divergence)
d_sqrt_jaccard_binary <- sqrt(d_jaccard_binary)

results <- do.call(rbind, list(
  evaluate_distance("Bray-Curtis", d_bray, md),
  evaluate_distance("Quantitative Jaccard (Ruzicka)", d_jaccard_quantitative, md),
  evaluate_distance("Binary Jaccard", d_jaccard_binary, md),
  evaluate_distance("phyloseq jsd (raw divergence)", d_js_divergence, md),
  evaluate_distance("sqrt(Jensen-Shannon divergence)", d_sqrt_js, md),
  evaluate_distance("sqrt(binary Jaccard)", d_sqrt_jaccard_binary, md)
))

print(results, digits = 8, row.names = FALSE)
