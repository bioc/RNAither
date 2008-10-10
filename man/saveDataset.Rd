\name{saveDataset}
\alias{saveDataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Save the normalized dataset into a dataset text file }
\description{
  Saves the normalized dataset and corresponding header into the specified dataset text file.
}
\usage{
saveDataset(header, data, dataSetFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{data}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{dataSetFile}{ the name of the text file the data will be saved in; can be the same as the old file (will be overwritten without prompting) }
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- divNorm(header, dataset, list(median, 2, 1, "SigIntensity", 1))
newheader <- normres[[1]]
newdataset <- normres[[2]]
saveDataset(newheader, newdataset, "save_testfile1.txt")
}

\keyword{ IO }
