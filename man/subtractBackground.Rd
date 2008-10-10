\name{subtractBackground}
\alias{subtractBackground}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Background substraction }
\description{
  Substracts a specified background value from the intensity values.
}
\usage{
subtractBackground(header, dataset, listOfArgs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing: 

- a character string specifying the column whose values will be used for background substraction

- a character string specifying the column whose values will be used as background
}}

\value{
  A list containing:
  \item{header }{The new header (with an added entry about the normalization procedure in the comments) }
  \item{dataset }{The new dataset with normalized values. The old values are saved in an extra column of the dataset with the suffix ".old" }
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- subtractBackground(header, dataset, list("SigIntensity", "Background"))
newheader <- normres[[1]]
newdataset <- normres[[2]]
}

\keyword{ manip }

