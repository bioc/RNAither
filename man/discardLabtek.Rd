\name{discardLabtek}
\alias{discardLabtek}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remove a complete plate from the analysis }
\description{
  Removes a plate/LabTek from the analysis by setting its spot type in the dataset to -1.
}
\usage{
discardLabtek(data, screenNr, labtekNr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{screenNr}{ the number of the experiment that contains the plate to discard }
  \item{labtekNr}{ the number of the plate to discard }
}

\value{
A new dataset that stil contains the specified plate/LabTek, but excludes it from the further analysis by setting its \code{SpotType}s to -1.
}

\seealso{ \code{\link{discardWells}} }
\examples{
data(exampleDataset, package="RNAither")

newdataset <- discardLabtek(dataset, 2, 2)
}

\keyword{ arith }