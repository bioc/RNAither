\name{BScore}
\alias{BScore}

\title{ BScore normalization }
\description{
  Normalization with BScores (see References).
}
\usage{
BScore(header, dataset, listOfArgs)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing:
  
  - a character string specifying the column whose values will be used for normalization 

  - a flag specifying whether controls should be excluded for the computation of the median polish (1) or not (0)
 }
}

\value{
  A list containing:
  \item{header }{The new header (with an added entry about the normalization procedure in the comments) }
  \item{dataset }{The new dataset with normalized values. The old values are saved in an extra column of the dataset with the suffix ".old" }
}
\references{ C. Brideau, B. Gunter, B. Pikounis, and A. Liaw. Improved statistical methods for hit selection in high-throughput screening. J Biomol Screen, 8:634-647, 2003 }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")
normres <- BScore(header, dataset, list("SigIntensity", 0))
newheader <- normres[[1]]
newdataset <- normres[[2]]
}

\keyword{manip}
