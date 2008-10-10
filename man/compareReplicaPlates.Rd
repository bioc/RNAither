\name{compareReplicaPlates}
\alias{compareReplicaPlates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare replica plates }
\description{
  Generates plots comparing the same plates in different experiments pairwise.
}
\usage{
compareReplicaPlates(header, dataset, plotTitle, col4val)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{plotTitle}{ the plot title }
  \item{col4val}{ a character string specifying the column whose values will be used for the plot }
}
\value{
For each plate, plots of pairwise comparisons between replicate intensities are generated and saved as a pdf file named after the experiment name specified in the header concatenated with the \code{ plotTitle }.
}

\seealso{ \code{\link{compareReplicates}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

compareReplicaPlates(header, dataset, "Comparison of replica plate", "SigIntensity")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ arith }
