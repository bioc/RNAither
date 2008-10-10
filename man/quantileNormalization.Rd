\name{quantileNormalization}
\alias{quantileNormalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quantile normalization }
\description{
  Quantile normalization (see References)
}
\usage{
quantileNormalization(header, dataset, listOfArgs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing: 

- 1 or 2, 1 meaning a normalization per experiment, 2 meaning a normalization per plate

- a character string specifying the column whose values will be used for normalization }
}

\value{
Returns a list, containing:

  \item{header }{the new header (with an added entry about the normalization procedure in the comments)}
  \item{dataset }{the new dataset with normalized values. The old values are saved in an extra column with the suffix ".old"}
}
\references{ B.M. Bolstad, R.A. Irizarry, M. Astrand, and T.P. Speed. A Comparison of Normalization Methods for High Density Oligonucleotide Array Data Based on Variance and Bias. Bioinformatics, 19(2): 185-193, 2003 }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- quantileNormalization(header, dataset, list(2, "SigIntensity"))
newheader <- normres[[1]]
newdataset <- normres[[2]]
}

\keyword{ manip }
