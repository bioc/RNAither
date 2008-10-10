\name{plotHistoPerscreen}
\alias{plotHistoPerscreen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a histogram of the data values per experiment }
\description{
  Plots and saves a histogram of the chosen data values.
}
\usage{
plotHistoPerscreen(header, dataset, channel, plotTitle, plotDesign, showPlot)
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
  \item{minOfScreens }{the number of the first experiment}
  \item{numOfScreens }{the number of the last experiment}
}

\seealso{ \code{\link{plotHisto}}, \code{\link{plotHistoPerplate}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

plotHistoPerscreen(header, dataset, "SigIntensity", "Distribution of Data and Controls", 1, 1)
}
\keyword{ device }
