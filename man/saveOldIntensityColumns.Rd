\name{saveOldIntensityColumns}
\alias{saveOldIntensityColumns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Save old intensity value columns }
\description{
  Duplicates the specified column and adds it to the end of the dataset.
}
\usage{
saveOldIntensityColumns(dataset, col4val)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{col4val}{ a character string specifying the column whose values will be saved as an extra column before normalization }
}

\value{
The values in the chosen column are saved in an extra column with the suffix ".old".
}

\examples{
data(exampleDataset, package="RNAither")
newdataset <- saveOldIntensityColumns(dataset, "SigIntensity")
}

\keyword{ datagen }
\keyword{ manip }
