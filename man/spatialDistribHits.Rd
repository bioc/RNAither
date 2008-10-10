\name{spatialDistribHits}
\alias{spatialDistribHits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting the spatial distribution of the hits }
\description{
  Plots the plates showing the spatial distribution of the hits using the \code{ plotPlate} function of the \code{ prada} package.
}
\usage{
spatialDistribHits(header, dataset, plotTitle, col4hits, col4anno, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{plotTitle}{ the plot title }
  \item{col4hits}{ a character vector specifying the name of the dataset column containing the binary hit vector }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"} }
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

\examples{
data(exampleHeader, package="RNAither")
data(pValVec1, package="RNAither")
data(scoredDataset1, package="RNAither")

##for details on the generation of pValVec1 and scoredDataset1, see \code{\link{Ttest}}

scoredHits1 <- hitselectionPval(scoredDataset1, pValVec1, "SigIntensity", "Hits1", 0.05, 
"GeneName", "pvalue_testfile1.txt")

hitDataset1 <- scoredHits1[[1]]

spatialDistribHits(header, hitDataset1, "Spatial distribution of hits", "Hits1", "GeneName", 1)
}

\keyword{ device }


