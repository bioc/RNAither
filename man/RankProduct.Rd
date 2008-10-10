\name{RankProduct}
\alias{RankProduct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform a Rank Product test }
\description{
  Performs the non-parametric rank product test on the intensity data.
}
\usage{
RankProduct(dataset, listofargs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listofargs}{ a list containing: 

- the number of permutations to perform to compute the p-values (usually 100)

- 1 or 2, depending if the search is for a significant decrease or increase

- a character string specifying the column whose values will be used for the test

- a character string specifying the name of the dataset column to be used to define the replicate, for example \code{ "GeneName"} or \code{ "Internal_GeneID"}
}}

\value{
Returns a list containing

  \item{pValVec }{a named vector of p-values}
  \item{dataset }{the dataset with an added column \code{ "p.value.rankproduct"}}
  \item{paste("pValue.rankproduct", testType, sep="_") }{the character string \code{"p.value.rankproduct"}}
  \item{"Rank product test" }{the character string "Rank product test"}
  
  The p values returned are equivalent to the percentage of false prediction (pfp), which in theory is the equivalent of false discovery rate (FDR). It is possible that they are larger than 1.
}

\seealso{ \code{\link{Ttest}}, \code{\link{MannWhitney}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

pvals1 <- RankProduct(dataset, list(100, 1, "SigIntensity", "GeneName"))
pValVec1 <- pvals1[[1]]
scoredDataset1 <- pvals1[[2]]
}

\keyword{ univar }
