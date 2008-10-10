\name{compareReplicateSD}
\alias{compareReplicateSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot the standard deviation of replicates }
\description{
  In the same fashion as \code{\link{spatialDistrib}}, generates a plot of the standard deviation of replicate values.
}
\usage{
compareReplicateSD(header, dataset, plotTitle, colname4SD, col4anno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{plotTitle}{ the plot title }
  \item{colname4SD}{ a character string specifying the column whose values will be used for the computation of the replicate standard deviation }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"} }
}

\value{
Generates a plot of the standard deviation of replicate values of all experiments. The plot is saved as a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

Wells showing positive controls sd are marked with a "P", wells showing negative controls sd with an "N".

The plot will also be saved as an html file containing mouse-overs with the siRNA name for each well.

The function returns the plotname.
}

\seealso{ \code{\link{spatialDistrib}}, \code{\link{compareReplicateSDPerScreen}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

compareReplicateSD(header, dataset, "Replicate standard intensity deviation", 
"SigIntensity", "GeneName")
}
\keyword{ arith }
\keyword{ device }