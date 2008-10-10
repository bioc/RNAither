\name{lowessNorm}
\alias{lowessNorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lowess normalization }
\description{
  Performs a plate-wise lowess normalization of the data.
}
\usage{
lowessNorm(header, dataset, listOfArgs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing:
  
- a character string specifying the column used as channel 1 (colname4ch1)

- a character string specifying the column used as channel 2 (colname4ch2)
}
}
\value{
Corrects intensity values in case the values of ch2 decrease with the increase of ch1 values.

Returns a list containing:
  \item{header }{the new header (with an added entry about the normalization procedure in the comments)}
  \item{dataset }{the new dataset with normalized values. The old values are saved in an extra column with the suffix ".old"}
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- lowessNorm(header, dataset, list("NbCells","SigIntensity"))
newheader <- normres[[1]]
newdataset <- normres[[2]]
}

\keyword{ manip }
