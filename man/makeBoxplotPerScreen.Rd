\name{makeBoxplotPerScreen}
\alias{makeBoxplotPerScreen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate a boxplot of the data per experiment }
\description{
  Generates a boxplot of intensity values per experiment.
}
\usage{
makeBoxplotPerScreen(header, dataset, channel, plotTitle, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the column whose values will be used for the boxplot }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
A boxplot of intensity values per experiment will be saved as a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.
}

\seealso{ \code{\link{makeBoxplotControls}}, \code{\link{makeBoxplotControlsPerPlate}}, \code{\link{makeBoxplotControlsPerScreen}}, \code{\link{makeBoxplotPerPlate}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

makeBoxplotPerScreen(header, dataset, "SigIntensity", "Data vs. Controls", 1)
}
\keyword{ device }
