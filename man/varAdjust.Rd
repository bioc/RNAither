\name{varAdjust}
\alias{varAdjust}

\title{ Variance adjustment }
\description{
  Divides the intensity values by their median absolute deviation (of the experiment or of the plate)
}
\usage{
varAdjust(header, dataset, listOfArgs)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing: 

- a character string specifying the column whose values will be used for normalization

- 1 or 2, 1 meaning a normalization per screen, 2 a normalization per plate

- a flag specifying whether controls should be excluded for the computation of the median absolute deviation (1) or not (0). }
}

\value{
  Divides the intensity values by their median absolute deviation (of the experiment or of the plate).
  
  Returns a list containing:
  \item{header }{The new header (with an added entry about the normalization procedure in the comments) }
  \item{dataset }{The new dataset with normalized values. The old values are saved in an extra column of the dataset with the suffix ".old" }
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- varAdjust(header, dataset, list(1, "SigIntensity", 1))
newheader <- normres[[1]]
newdataset <- normres[[2]]
}

\keyword{ manip }
