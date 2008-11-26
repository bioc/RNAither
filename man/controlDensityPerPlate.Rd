\name{controlDensityPerPlate}
\alias{controlDensityPerPlate}

\title{ Plotting the control density per plate }
\description{
  Plots the density of positive and negative controls (if applicable) for each plate.
}
\usage{
controlDensityPerPlate(header, dataset, channel, plotTitle, plotDesign, showPlot, supHisto)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the name of the column containing the values for computing the density, e.g. \code{"SigIntensity"} }
  \item{plotTitle}{ the plot title }
  \item{plotDesign}{ 1 or 2. 1 will generate one window containing all plots, 2 will generate a series of plots }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
  \item{supHisto}{ 0 or 1. 1 will additionally superimpose a colour histogram of the values for the positive and negative controls. Otherwise choose 0.}
}

\value{

Generates a series of plots for each experiment and each plate, showing the density of positive and negative controls (if applicable). The plots are saved as pdf and png files named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns a list of length 3 containing: 

  \item{plotName }{ the plot name }
  \item{c(minOfScreens, numOfScreens) }{ a vector with the number of the first experiment and of the last experiment }
  \item{c(minOfPlates, numOfPlates) }{ a vector with the number of the first plate and the number of the last plate }
}

\seealso{ \code{\link{controlDensity}}, \code{\link{controlDensityPerScreen}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

controlDensityPerPlate(header, dataset, "SigIntensity", "Control density", 1, 1, 1)
}
\keyword{ device }
