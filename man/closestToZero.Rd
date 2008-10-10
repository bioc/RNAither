\name{closestToZero}
\alias{closestToZero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Return the replicate value closest to zero }
\description{
  Out of a set of replicate values, returns the one closest to zero.
}
\usage{
closestToZero(Ivec, na.rm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ivec}{ All channel values for a specific siRNA/gene }
  \item{na.rm}{ Removes NA values }
}

\value{
  A double giving the value closest to zero out of the given replicate values.
}

\seealso{ \code{\link{rms}}, \code{\link{trim}}, \code{\link{furthestFromZero}}, \code{\link{summarizeReps}}, \code{\link{summarizeRepsNoFiltering}} }

\examples{
data(exampleDataset, package="RNAither")
Indexes <- findReplicates(dataset, "GeneName", "CPSF1")
replicateclosest <- closestToZero(dataset$SigIntensity[Indexes])
}
\keyword{arith}
