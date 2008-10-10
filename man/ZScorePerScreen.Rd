\name{ZScorePerScreen}
\alias{ZScorePerScreen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ZScore normalization per experiment }
\description{
  ZScore normalization not per plate, but per experiment (see Value and References)
}
\usage{
ZScorePerScreen(header, dataset, listOfArgs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing:

- a character string specifying the column whose values will be used for normalization

- a flag specifying whether controls should be excluded for the computation of the median and median absolute deviation (1) or not (0). }
}
\value{
 The ZScore is defined as the quotient of the difference between an intensity value and the median of the experiment, and of the median absolute deviation.
 
 Returns a list containing:
   \item{header }{The new header (with an added entry about the normalization procedure in the comments) }
   \item{dataset }{The new dataset with normalized values. The old values are saved in an extra column of the dataset with the suffix ".old" }
   }
   
\references{ N. Malo et al. Statistical practice in high-throughput screening data analysis. Nature Biotech, 24(2): 167-175, 2006. }


\seealso{ \code{\link{ZScore}}, \code{\link{BScore}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- ZScorePerScreen(header, dataset, list("SigIntensity", 1))
newheader <- normres[[1]]
newdataset <- normres[[2]]
}

\keyword{ manip }