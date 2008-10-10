\name{trim}
\alias{trim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the replicate mean with trimmed values }
\description{
  Computes the mean of replicate values, omitting the highest and the lowest 5%.
}
\usage{
trim(Ivec, na.rm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ivec}{ All channel values for a specific siRNA/gene }
  \item{na.rm}{ Removes NA values }
}

\value{
A double giving the trimmed mean of the given replicate values, i.e. omitting the highest and the lowest 5%.
}

\seealso{ \code{\link{rms}}, \code{\link{closestToZero}}, \code{\link{furthestFromZero}}, \code{\link{summarizeReps}}, \code{\link{summarizeRepsNoFiltering}} }


\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

Indexes <- findReplicates(dataset, "GeneName", "CPSF1")
replicatemean <- trim(dataset$SigIntensity[Indexes])
}

\keyword{ arith }
