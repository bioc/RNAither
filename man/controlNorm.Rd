\name{controlNorm}
\alias{controlNorm}

\title{ Normalization on controls }
\description{
  Performs a normalization on either positive or negative controls.
}
\usage{
controlNorm(header, dataset, listOfArgs)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing: 
  
    - a character string specifying the column whose values will be used for normalization
    
    - 1 or 2, 1 meaning a normalization per experiment, 2 meaning a normalization per plate
    
    - 0 or 1, 0 meaning a normalization on the median of negative controls, 1 meaning a normalization on the median of positive controls. Can also be the \code{ GeneName} of a specific control siRNA
    
    - 1 or 2, 1 meaning the signal values are divided by the median, 2 meaning the median is subtracted from the signal values
  }
}

\value{

Returns a list containing:

  \item{header }{ the new header (with an added entry about the normalization procedure in the comments). }
  \item{dataset }{ the new dataset with normalized values. The old values are saved in an extra column with the suffix ".old". }
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- controlNorm(header, dataset, list(2, 0, "SigIntensity", 1))

newheader <- normres[[1]]
newdataset <- normres[[2]]
}

\keyword{ manip }

