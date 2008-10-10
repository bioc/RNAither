\name{compareReplicates}
\alias{compareReplicates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare replicate values }
\description{
  Plots replicate intensities pairwise for each experiment.
}
\usage{
compareReplicates(header, dataset, plotTitle, col4val, col4anno, plotDesign, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{plotTitle}{ the plot title }
  \item{col4val}{ a character string specifying the column whose values will be used for the plot }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"} }
  \item{plotDesign}{ 1 or 2. 1 will generate one window containing all plots, 2 will generate a series of plots }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
}

\value{
For each experiment, plots of pairwise comparisons between replicate intensities are generated and saved as a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}, and the number of the experiment.

The function returns a list containing:
  \item{plotName }{ the plot name }
  \item{minOfScreens }{ the number of the first experiment }
  \item{numOfScreens }{ the number of the last experiment }
  \item{maxCombinationNum }{ the number of replicates to compare }
}

\seealso{ \code{\link{compareReplicaPlates}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

compareReplicates(header, dataset, "Comparison of Replicates", "SigIntensity", "GeneName", 1, 0)
}

\keyword{ arith }
\keyword{ device }

