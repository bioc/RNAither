\name{controlDensity}
\alias{controlDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting the control density }
\description{
  Plots the density of positive and negative controls (if applicable) for all controls contained in the dataset.
}
\usage{
controlDensity(header, dataset, channel, plotTitle, showPlot)
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
Plots the density of positive and negative controls (if applicable) for all controls contained in the dataset. The plot is saved as a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns the plot name.
}

\seealso{ \code{\link{controlDensityPerScreen}}, \code{\link{controlDensityPerPlate}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

controlDensity(header, dataset, "SigIntensity", "Control density", 1)
}

\keyword{ device }
