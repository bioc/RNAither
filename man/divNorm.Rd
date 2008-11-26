\name{divNorm}
\alias{divNorm}

\title{ Mean, median, ... , normalization }
\description{
  Normalization with the mean, median, or any other function.
}
\usage{
divNorm(header, dataset, listOfArgs)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing: 
  
- a character string specifying the column whose values will be used for normalization

- a function to be used for the normalization, e.g. \code{ mean}, \code{ median}, ...

- 1 or 2, 1 meaning a normalization per experiment, 2 meaning a normalization per plate

- 1 or 2, 1 meaning the normalization is achieved by a division of the intensity values by the outcome of \code{ funname}, 2, meaning by a substraction

- a flag specifying whether controls should be excluded for the computation of the result of the function specified in the first element (1) or not (0).
  }
}

\value{
Returns a list containing: 

  \item{header }{ the new header (with an added entry about the normalization procedure in the comments) }
  \item{dataset }{ the new dataset with normalized values. The old values are saved in an extra column with the suffix ".old" }
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- divNorm(header, dataset, list(median, 2, 1, "SigIntensity", 1))

newheader <- normres[[1]]
newdataset <- normres[[2]]

}

\keyword{ manip }
