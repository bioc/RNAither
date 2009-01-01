\name{compareHits}
\alias{compareHits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Searching for common hits between different scoring methods }
\description{
  Searches for common hits between different scoring methods.
}
\usage{
compareHits(hitVec1, hitVec2, namesHitVec1, namesHitVec2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hitVec1, hitVec2}{ the two binary hit vectors to be compared }
  \item{namesHitVec1, namesHitVec2}{ the names of the siRNAs corresponding to the hit vectors }
}
\value{
Returns a character vector indicating which siRNAs are identified as hits in two different hit scoring schemes.
}

\seealso{ \code{\link{vennDiag}} , \code{\link{Ttest}}, \code{\link{MannWhitney}} }
\examples{
data(scoredDataset1, package="RNAither")
data(pValVec1, package="RNAither")

data(scoredDataset2, package="RNAither")
data(pValVec2, package="RNAither")

##for details on the generation of pValVec and scoredDataset,
##see the examples of the functions Ttest and MannWhitney linked above.

scoredHits1 <- hitselectionPval(scoredDataset1, pValVec1, "SigIntensity", "Hits1", 0.05, 
"GeneName", "pvalue_testfile1.txt")
scoredHits2 <- hitselectionPval(scoredDataset2, pValVec2, "SigIntensity", "Hits2", 0.05, 
"GeneName", "pvalue_testfile2.txt")

hitVector1 <- scoredHits1[[2]]
hitVector2 <- scoredHits2[[2]]

common_hits <- compareHits(hitVector1, hitVector2, names(hitVector1), names(hitVector2))
}

\keyword{ arith }
