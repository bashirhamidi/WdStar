#' Calculate the Tw2 Statistic for Heteroscedastic Test
#'
#' This function calculates the \eqn{T_w^2} statistic to compare two groups
#' with potentially unequal sample sizes. It is suitable for microbiome and
#' other multivariate distance data.
#'
#' @param dm A distance matrix, representing dissimilarity between observations.
#' @param f A factor variable indicating the group for each observation.
#'
#' @return The calculated Tw2 statistic.
#'
#' @details
#' The function first checks if the factor variable has exactly two levels.
#' The \eqn{T_w^2} statistic is a modification of Hotelling's T-square statistic
#' adapted for heteroscedasticity in multivariate distance data. This function
#' returns the observed statistic only; use \code{\link{Tw2.test}} for
#' permutation-based significance testing.
#'
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
#'   Tw2(dm_12, f_12)
#' }
#'
#' @references
#' Alekseyenko AV. Multivariate Welch t-test on distances. \emph{Bioinformatics}.
#' 2016;32(23):3552-3558. doi:10.1093/bioinformatics/btw524
#'
#' @seealso \code{\link{Tw2.test}}
#'
#' @export
Tw2 <- function(dm, f) {
  if (nlevels(f) != 2) {
    return(NULL)
  }
  SS2 <- dist.ss2(as.matrix(dm^2), f)
  ns <- summary(f)
  N <- sum(ns)

  SST <- sum(SS2) / N
  SSW1 <- SS2[1, 1] / ns[1]
  SSW2 <- SS2[2, 2] / ns[2]
  SSW <- SSW1 + SSW2

  s1 <- SSW1 / (ns[1] - 1)
  s2 <- SSW2 / (ns[2] - 1)

  t.stat <- (ns[1] + ns[2]) / (ns[1] * ns[2]) * (SST - SSW) / (s1 / ns[1] + s2 / ns[2])
  t.stat
}

#' Calculate the Wd* Statistic for Heteroscedastic Test
#'
#' This function calculates the \eqn{\mathnormal{W}_d^*} statistic to compare
#' k groups with potentially unequal multivariate dispersions and sample sizes.
#' It is suitable for microbiome and other multivariate distance data.
#'
#' @param dm A distance matrix, representing dissimilarity between observations.
#' @param f A factor variable indicating the group for each observation.
#'
#' @return The calculated Wd* statistic.
#'
#' @details
#' This function is an extension of Welch's ANOVA statistic suitable for
#' multivariate distance data. The \eqn{\mathnormal{W}_d^*} statistic is
#' computed from pairwise squared distances among groups while accounting for
#' unequal group sizes and differences in multivariate spread across groups.
#' This function returns the observed statistic only; use
#' \code{\link{WdS.test}} for permutation-based significance testing.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("enterotype", package = "phyloseq")
#'   ent <- phyloseq::subset_samples(enterotype, !is.na(Enterotype))
#'   dm <- phyloseq::distance(ent, method = "bray")
#'   meta <- data.frame(phyloseq::sample_data(ent))
#'   f <- factor(meta$Enterotype)
#'
#'   WdS(dm, f)
#' }
#'
#' @references
#' Hamidi, Bashir, et al. "$ W_ {d}^{*} $-test: robust distance-based
#' multivariate analysis of variance." Microbiome 7.1 (2019): 1-9.
#'
#' @seealso \code{\link{WdS.test}}, \code{\link{Tw2}},
#'   \url{https://github.com/alekseyenko/WdStar}
#'
#' @export
WdS <- function(dm, f) {
  ns <- table(f)
  SS2 <- dist.ss2(as.matrix(dm)^2, f)
  s2 <- diag(SS2) / ns / (ns - 1)
  W <- sum(ns / s2)

  idxs <- apply(utils::combn(levels(f), 2), 2, function(idx) levels(f) %in% idx)

  Ws <- sum(apply(
    idxs, 2,
    function(idx) {
      sum(ns[idx]) / prod(s2[idx]) *
        (sum(SS2[idx, idx]) / sum(ns[idx]) - sum(diag(SS2[idx, idx]) / ns[idx]))
    }
  ))
  k <- nlevels(f)
  h <- sum((1 - ns / s2 / W)^2 / (ns - 1))
  Ws / W / (k - 1) / (1 + (2 * (k - 2) / (k^2 - 1)) * h)
}
