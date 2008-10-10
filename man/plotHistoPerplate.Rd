\name{plotHistoPerplate}
\alias{plotHistoPerplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a histogram of the data values per plate }
\description{
  Plots and saves a histogram of the chosen data values per experiment and per plate.
}
\usage{
plotHistoPerplate(header, dataset, channel, plotTitle, plotDesign, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the name of the column containing the values to be plotted, e.g. \code{ "SigIntensity"} }
  \item{plotTitle}{ the plot title }
  \item{plotDesign}{ 1 or 2. 1 will generate one window containing all plots, 2 will generate a series of plots. }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
Saves the histograms in a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns a list containing:

  \item{histoName }{the plotname}
  \item{c(minOfScreens, numOfScreens) }{a vector with the number of the first experiment and of the last experiment}
  \item{c(minOfPlates, numOfPlates) }{a vector with the number of the first plate and the number of the last plate}
}

\seealso{ \code{\link{plotHisto}}, \code{\link{plotHistoPerscreen}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

plotHistoPerplate(header, dataset, "SigIntensity", "Distribution of Data and Controls", 1, 1)
}
\keyword{ device }
