\name{MannWhitney}
\alias{MannWhitney}

\title{ Perform a Mann-Whitney test }
\description{
  Performs the non-parametric Mann-Whitney test on the intensity data.
}
\usage{
MannWhitney(dataset, listofargs)
}

\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listofargs}{ a list containing: 

- \code{"g"} (greater) for significant increase, \code{"l"} (lower) for significant decrease, or \code{"two.sided"} for both 

- either a number indicating the true value of the mean, or a character string indicating the name of the gene to compare with 

- a character string specifying the column whose values will be used for the test

- a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"}
}
}

\value{
Returns a list containing:

  \item{pValVec }{a named vector of p-values}
  \item{dataset }{the dataset with an added column \code{ "p.value.mannwhitney"}}
  \item{paste(\code{"pValue.mannwhitney"}, testType, sep="_") }{the character string "p.value.mannwhitney" concatenated with the \code{testType} (first element of \code{listofargs})}
  \item{"Mann-Whitney test" }{the character string "Mann-Whitney test"}
}

\seealso{ \code{\link{Ttest}}, \code{\link{RankProduct}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

pvals1 <- MannWhitney(dataset, list("l", median(dataset$SigIntensity, na.rm=TRUE), 
"SigIntensity", "GeneName"))
pValVec1 <- pvals1[[1]]
scoredDataset1 <- pvals1[[2]]
}

\keyword{ manip }

