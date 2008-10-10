\name{eraseDataSetColumn}
\alias{eraseDataSetColumn}

\title{ Remove columns from dataset }
\description{
  Removes a specified column from a dataset.
}
\usage{
eraseDataSetColumn(dataset, colname)
}

\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{colname}{ a character string specifying the name of the column to be removed }
}

\value{
  An R data frame with the specified column removed.
}

\examples{
data(exampleDataset, package="RNAither")

newdataset <- eraseDataSetColumn(dataset, "SDSIntensity")
}

\keyword{ manip }

