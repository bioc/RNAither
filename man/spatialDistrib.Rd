\name{spatialDistrib}
\alias{spatialDistrib}

\title{ Generate spatial plots of intensity values }
\description{
  Generate plots of plates and their intensity values.
}
\usage{
spatialDistrib(header, dataset, plotTitle, col4plot, col4anno, showPlot)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{plotTitle}{ the plot title }
  \item{col4plot}{ a character string specifying the column whose values will be used for the plot }
  \item{col4anno}{ a character string specifying the column whose values will be used for the annotation of the plot }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
For each plate, the plot will be saved as a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}, the number of the experiment, and the number of the plate.

Wells containing positive controls are marked with a "P", wells containing negative controls with an "N".

Each plate will also be saved as an html file containing mouse-overs with the siRNA name for each well.


The function returns a list containing:

  \item{histoName }{the plotname}
  \item{c(minOfScreens, numOfScreens) }{a vector with the number of the first experiment and of the last experiment}
  \item{c(minOfPlates, numOfPlates) }{a vector with the number of the first plate and the number of the last plate}
}

\seealso{ \code{\link{compareReplicateSD}}, \code{\link{compareReplicateSDPerScreen}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

spatialDistrib(header, dataset, "Spatial distribution of cell counts", "NbCells", "GeneName", 1)
}

\keyword{ device }

