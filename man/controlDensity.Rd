\name{controlDensity}
\alias{controlDensity}

\title{ Plotting the control density }
\description{
  Plots the density of positive and negative controls (if applicable) for all controls contained in the dataset.
}
\usage{
controlDensity(header, dataset, channel, plotTitle, showPlot, supHisto)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the name of the column containing the values for computing the density, e.g. \code{"SigIntensity"} }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
  \item{supHisto}{ 0 or 1. 1 will additionally superimpose a colour histogram of the values for the positive and negative controls. Otherwise choose 0.}
}

\value{
Plots the density of positive and negative controls (if applicable) for all controls contained in the dataset. The plot is saved as a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns the plot name.
}

\seealso{ \code{\link{controlDensityPerScreen}}, \code{\link{controlDensityPerPlate}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

controlDensity(header, dataset, "SigIntensity", "Control density", 1, 1)
}

\keyword{ device }
