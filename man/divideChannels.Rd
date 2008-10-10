\name{divideChannels}
\alias{divideChannels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Divide channel values }
\description{
  Replace two channels by their ratio.
}
\usage{
divideChannels(ch1, ch2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ch1}{ a vector giving all values from channel 1 }
  \item{ch2}{ a vector giving all values from channel 2 }
}

\value{
A vector of the ratio of channel 1 and channel 2.
}

\seealso{ \code{\link{sumChannels}} }
\examples{
data(exampleDataset, package="RNAither")

newch <- divideChannels(dataset$SigIntensity, dataset$NbCells)
}

\keyword{ manip }
\keyword{ arith }
