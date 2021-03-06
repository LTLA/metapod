% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groupedFisher.R
\name{groupedFisher}
\alias{groupedFisher}
\title{Combine grouped p-values with Fisher's method}
\usage{
groupedFisher(p.values, grouping, log.p = FALSE)
}
\arguments{
\item{p.values}{A numeric vector containing p-values for individual tests.}

\item{grouping}{A vector or factor of length equal to \code{p.values}, specifying the group to which each test is assigned.

Alternatively, an \link{rle} object where each run corresponds to a group and specifies the entries of \code{p.values} belonging to that group.
This assumes that \code{p.values} is ordered such that all entries in the same group are adjacent to each other.}

\item{log.p}{Logical scalar indicating whether the p-values in \code{p.values} are log-transformed.}
}
\value{
A list containing:
\itemize{
\item \code{p.value}, a named numeric vector of length equal to the number of unique levels in \code{grouping}.
This contains the Fisher combined p-value for each group, log-transformed if \code{log.p=TRUE}.
Each entry is named according to the group.
\item \code{representative}, a named integer scalar specifying the representative test for each group.
Each index refers to an entry of \code{p.values} and is named according to its group.
\item \code{influential}, a logical vector of length equal to \code{p.values}.
Entries are \code{TRUE} for any p-value that is deemed \dQuote{influential} to the final combined p-value for its group.
}
}
\description{
Combine p-values from grouped tests with Fisher's method.
Groups are defined according to unique levels of a grouping factor.
}
\details{
The joint null hypothesis for each group is that all of the individual null hypotheses are true.
Fisher's method combines information from all individual nulls to determine if the joint null should be rejected.
Compared to Stouffer's and Pearson's methods, Fisher's method provides more sensitivity to the smallest individual p-value.
This method is only applicable to independent tests and no weights are considered.

The representative test for each group is defined as the test with the lowest p-value, as this has the greatest effect on the combined p-value. 
All tests for each group are considered to be influential as increasing any of them (e.g., to unity) would result in a larger combined p-value.
}
\examples{
p1 <- rbeta(100, 0.8, 1)
g <- sample(10, length(p1), replace=TRUE)

# Standard application:
out <- groupedFisher(p1, g)
str(out)

# With log p-values. 
out <- groupedFisher(log(p1), g, log.p=TRUE)
str(out)

}
\references{
Fisher RA (1925).
\emph{Statistical Methods for Research Workers}.
Oliver and Boyd (Edinburgh).
}
\seealso{
\code{\link{parallelFisher}}, for a version that operates on parallel vectors of p-values.

\code{\link{groupedStouffer}} and \code{\link{groupedPearson}}, for different approaches to testing a joint null of independent hypotheses.
}
\author{
Aaron Lun
}
