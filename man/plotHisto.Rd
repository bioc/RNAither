\name{plotHisto}
\alias{plotHisto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a histogram of the data values }
\description{
  Plots and saves a histogram of the chosen data values.
}
\usage{
plotHisto(header, dataset, channel, plotTitle, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the name of the column containing the values to be plotted, e.g. \code{ "SigIntensity"} }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
Saves the histogram in a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns the plot name.
}

\seealso{ \code{\link{plotHistoPerplate}}, \code{\link{plotHistoPerscreen}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

plotHisto(header, dataset, "SigIntensity", "Distribution of Data and Controls", 1)
}
\keyword{ device }
