\name{joinDatasets}
\alias{joinDatasets}

\title{ Join datasets }
\description{
  Merges two or more datasets into one.
}
\usage{
joinDatasets(listOfDatasets)
}

\arguments{
  \item{listOfDatasets}{ a list of the datasets to join }
}
\value{
  The joined datasets.
}

\seealso{ \code{\link{generateDatasetFile}}, \code{\link{joinDatasetFiles}} }

\examples{
data(exampleDataset, package="RNAither")
doubledataset <- joinDatasets(list(dataset, dataset))

}

\keyword{ datagen }

