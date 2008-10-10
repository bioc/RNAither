\name{controlDensityPerScreen}
\alias{controlDensityPerScreen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting the control density per experiment }
\description{
  Plots the density of positive and negative controls (if applicable) for each experiment.
}
\usage{
controlDensityPerScreen(header, dataset, channel, plotTitle, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the name of the column containing the values for computing the density, e.g. \code{"SigIntensity"} }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
}

\value{
Generates a series of plots for each experiment, showing the density of positive and negative controls (if applicable). The plots are saved as pdf and png files named after the experiment name specified in the header concatenated with the \code{ plotTitle} and the number of the experiment.

The function returns a list of length 3 containing:

  \item{plotName }{ the plotname }
  \item{minOfScreens }{ the number of the first experiment }
  \item{numOfScreens }{ the number of the last experiment }
}

\seealso{ \code{\link{controlDensity}}, \code{\link{controlDensityPerPlate}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

controlDensityPerScreen(header, dataset, "SigIntensity", "Control density", 1)
}
\keyword{ device }
