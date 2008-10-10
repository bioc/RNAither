\name{discardWells}
\alias{discardWells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remove wells from the analysis }
\description{
  Removes wells from the analysis by setting their spot type in the dataset to -1.
}
\usage{
discardWells(data, screenNr, labtekNr, vecPositions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{screenNr}{ the number of the experiment that contains the plate to discard }
  \item{labtekNr}{ the number of the plate to discard }
  \item{vecPositions}{ a vector specifying the numbers of the wells to discard }
}

\value{
A new dataset that does not contain the specified wells.
A new dataset that stil contains the specified wells/spots, but excludes them from the further analysis by setting their \code{SpotType}s to -1.
}

\seealso{ \code{\link{discardLabtek}} }
\examples{
data(exampleDataset, package="RNAither")

newdataset <- discardWells(dataset, 2, 1, c(1, 10, 15))
}

\keyword{ arith }
