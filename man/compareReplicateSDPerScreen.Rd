\name{compareReplicateSDPerScreen}
\alias{compareReplicateSDPerScreen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot the standard deviation of replicates for each experiment}
\description{
  In the same fashion as \code{\link{spatialDistrib}}, generates plots of the standard deviation of replicate values for each experiment.
}
\usage{
compareReplicateSDPerScreen(header, dataset, plotTitle, colname4SD, col4anno, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{plotTitle}{ the plot title }
  \item{colname4SD}{ a character string specifying the column whose values will be used for the computation of the replicate standard deviation }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"} }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
}

\value{

Generates plots of the standard deviation of replicate values for each experiment. The plots are saved as png files named after the experiment name specified in the header concatenated with the \code{ plotTitle} and the number of the experiment.

Wells showing positive controls sd are marked with a "P", wells showing negative controls sd with an "N".

The plots will also be saved as html files containing mouse-overs with the siRNA name for each well.

The function returns a list of length 3 containing: 

  \item{basicPlotName }{the plot name}
  \item{minOfScreens }{the number of the first experiment}
  \item{numOfScreens }{the number of the last experiment}
  
}

\seealso{ \code{\link{spatialDistrib}}, \code{\link{compareReplicateSD}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

compareReplicateSDPerScreen(header, dataset, "Replicate standard intensity deviation", 
"SigIntensity", "GeneName", 1)
}

\keyword{ arith }
\keyword{ device }
