#' Combine parallel p-values with Simes' method
#'
#' Combine p-values from parallel tests with Simes' method.
#' Each group of p-values is defined from the corresponding entries across all vectors.
#'
#' @param p.values A list of numeric vectors of the same length, containing the p-values to be combined.
#' @param weights A numeric vector of positive weights, with one value per vector in \code{...}.
#' Each weight is applied to all entries of itscorresponding vector, i.e., all p-values in that vector receive the same weight.
#'
#' Alternatively, a list of numeric vectors of weights with the same structure as \code{p.values}.
#' Each p-value is then assigned the weight in the corresponding entry of \code{weights}.
#'
#' Alternatively \code{NULL}, in which case all p-values are assigned equal weight.
#' @param log.p Logical scalar indicating whether the p-values in \code{p.values} are log-transformed.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{p.value}, a numeric vector of length equal to the length of each vector in \code{p.values}.
#' This contains the combined Simes p-value for each group, log-transformed if \code{log.p=TRUE}.
#' \item \code{representative}, an integer scalar specifying the representative test in each group.
#' \item \code{influential}, a list of logical vectors mirroring the structure of \code{p.values}.
#' Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value.
#' }
#' 
#' @details
#' The joint null hypothesis for each group is that all of the individual null hypotheses are true.
#' Simes' method will reject the joint null if any of the individual nulls are rejected, providing weak control of the family-wise error rate.
#' 
#' In theory, the method is only applicable to independent tests, but experience suggests that it is quite robust to dependencies.
#' The calculation itself is very closely related to the Benjamini-Hochberg method for controlling the false discovery rate.
#' One can actually obtain Simes' combined p-value by taking the smallest BH-adjusted p-value across a group.
#'
#' If non-equal weights are provided, they are treated as relative frequency weights.
#' That is, if one p-value is given a weight of 10 and another p-value is given a weight of 1, 
#' the former is considered to occur 10 times more frequently than the latter.
#'
#' The representative test for each group is defined as the test with the p-value that is ultimately used as the combined p-value.
#' Briefly, one can identify this test as that with the smallest BH-adjusted p-value if the monotonicity adjustment were omitted.
#' The influential tests for each group are defined as those with p-values no greater than the representative test's p-value.
#' This is based on the fact that increasing them (e.g., by setting them to unity) would result in a larger combined p-value.
#' 
#' @author Aaron Lun
#' @examples
#' p1 <- rbeta(100, 0.8, 1)
#' p2 <- runif(100)
#' p3 <- rbeta(100, 0.5, 1)
#'
#' # Standard application:
#' out <- parallelSimes(list(p1, p2, p3))
#' str(out)
#'
#' # With vector-level weights:
#' out <- parallelSimes(list(p1, p2, p3), weights=c(10, 20, 30))
#' str(out)
#' 
#' # With log p-values. 
#' out <- parallelSimes(list(log(p1), log(p2), log(p3)), log.p=TRUE)
#' str(out)
#'
#' @references
#' Simes RJ (1986).
#' An improved Bonferroni procedure for multiple tests of significance.
#' \emph{Biometrika} 73:751-754.
#'
#' Sarkar SK and Chung CK (1997).
#' The Simes method for multiple hypothesis testing with positively dependent test statistics.
#' \emph{J. Am. Stat. Assoc.} 92, 1601-1608.
#' 
#' Benjamini Y and Hochberg Y (1997).
#' Multiple hypotheses testing with weights.
#' \emph{Scand. J. Stat.} 24, 407-418.
#'
#' @seealso
#' \code{\link{groupedSimes}}, for a version that combines p-values based on a grouping factor.
#' @export
parallelSimes <- function(p.values, weights=NULL, log.p=FALSE) {
    compute_parallel_simes(p.values, weights, log.p)
}
