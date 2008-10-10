\name{Ttest}
\alias{Ttest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform a Student's t-test }
\description{
  Performs a Student's t-test on the intensity data.
  }
\usage{
Ttest(dataset, listofargs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listofargs}{ a list containing: 

- \code{ "g"} (greater) for significant increase, \code{ "l"} (lower) for significant decrease, or \code{ "two.sided"} for both 

- either a number indicating the true value of the mean, or a character string indicating the name of the gene to compare with 

- a character string specifying the column whose values will be used for the test

- a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"}
}
}

\value{
Returns a list containing:

  \item{pValVec }{a named vector of p-values}
  \item{dataset }{the dataset with an added column \code{ "p.value.mannwhitney"}}
  \item{paste("pValue.ttest", testType, sep="_") }{the character string \code{"pValue.ttest"} concatenated with the \code{testType} (first element of \code{listofargs})}
  \item{"t test" }{the character string \code{"t test"}}
}

\seealso{ \code{\link{MannWhitney}}, \code{\link{RankProduct}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

pvals1 <- Ttest(dataset, list("l", median(dataset$SigIntensity, na.rm=TRUE), "SigIntensity", "GeneName"))
pValVec1 <- pvals1[[1]]
scoredDataset1 <- pvals1[[2]]
}

\keyword{ univar }