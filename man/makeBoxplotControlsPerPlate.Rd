\name{makeBoxplotControlsPerPlate}
\alias{makeBoxplotControlsPerPlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate a boxplot of the data vs. the controls for each plate }
\description{
  Generates a boxplot of intensity values of negative controls, positive controls and experimental data for each plate of each experiment available in the dataset.
}
\usage{
makeBoxplotControlsPerPlate(header, dataset, channel, plotTitle, plotDesign, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the column whose values will be used for the boxplot }
  \item{plotTitle}{ the plot title }
  \item{plotDesign}{ 1 or 2. 1 will generate one window containing all plots, 2 will generate a series of plots. }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
For each experiment, a series of boxplots of intensity values of negative controls, positive controls and experimental data will be saved as a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle} and the number of the experiment.

The function returns a list containing:

  \item{plotName }{the plotname}
  \item{c(minOfScreens, numOfScreens) }{a vector with the number of the first experiment and of the last experiment}
  \item{c(minOfPlates, numOfPlates) }{a vector with the number of the first plate and the number of the last plate}
}

\seealso{ \code{\link{makeBoxplotControls}}, \code{\link{makeBoxplotControlsPerScreen}}, \code{\link{makeBoxplotPerPlate}}, 
\code{\link{makeBoxplotPerScreen}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

makeBoxplotControlsPerPlate(header, dataset, "SigIntensity", "Data vs. Controls", 1, 1)
}

\keyword{ device }