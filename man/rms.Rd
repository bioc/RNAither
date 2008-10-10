\name{rms}
\alias{rms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the replicate root mean square }
\description{
  Computes the root mean square of replicate values
}
\usage{
rms(Ivec, na.rm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ivec}{ All channel values for a specific siRNA/gene }
  \item{na.rm}{ Removes NA values }
}

\value{
A double giving the root mean square of the given replicate values.
}

\seealso{ \code{\link{trim}}, \code{\link{closestToZero}}, \code{\link{furthestFromZero}}, \code{\link{summarizeReps}}, \code{\link{summarizeRepsNoFiltering}} }
\examples{
data(exampleDataset, package="RNAither")

Indexes <- findReplicates(dataset, "GeneName", "CPSF1")
rmsval <- rms(dataset$SigIntensity[Indexes])
}

\keyword{arith}