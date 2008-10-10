\name{channelPlot}
\alias{channelPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot signal channels against each other }
\description{
  Generates plots allowing pairwise comparison of signal channels. Fits a lowess regression curve into the plots.
}

\usage{
channelPlot(header, dataset, vecOfChannels, flag, plotTitle, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{vecOfChannels}{ A vector containing the names of the signal channels to be compared, e.g. "SigIntensity" }
  \item{flag}{ 0, 1, or 2. 0 uses the data from the complete dataset, 1 makes comparisons for each experiment, 2 makes comparisons for each plate. }
  \item{plotTitle}{ The plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
  Saves the plots in pdf and png files named after the experiment name specified in the header concatenated with the \code{plotTitle}, the number of the comparison, and if applicable the experiment number and/or the plate number.


  When \code{flag} == 0, returns the plot name (\code{plotName}). 

  
  When \code{flag} == 1, returns a list containing:
  \item{plotName }{ The plot name }
  \item{minOfScreens}{ The number of the first experiment }
  \item{numOfScreens}{ The number of the last experiment }
  When \code{flag} == 2, returns a list containing: the plot name, a vector with the number of the first experiment and of the last experiment, and a vector with the number of the first plate and the number of the last plate.
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")
plotname <- channelPlot(header, dataset, c("SigIntensity", "NbCells"), 0, "Channel comparison", 1)
}

\keyword{ device }