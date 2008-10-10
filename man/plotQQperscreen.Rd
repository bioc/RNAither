\name{plotQQperscreen}
\alias{plotQQperscreen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make a QQ plot per experiment }
\description{
  Shows and saves a QQ plot of the data for each experiment in the dataset.
}
\usage{
plotQQperscreen(header, dataset, channel, plotTitle, plotDesign, showPlot)
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
Saves the QQ plots in a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns a list containing: 

  \item{histoName }{the plotname}
  \item{minOfScreens }{the number of the first experiment}
  \item{numOfScreens }{the number of the last experiment}
}

\seealso{ \code{\link{plotQQ}}, \code{\link{plotQQperplate}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

plotQQperscreen(header, dataset, "SigIntensity", "QQplot", 1, 1)
}

\keyword{ device }
