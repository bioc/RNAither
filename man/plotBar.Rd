\name{plotBar}
\alias{plotBar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot signal intensities per well }
\description{
  Plots signal intensity values for each well, a blue line showing the median, two green lines showing one median absolute deviation, two red lines showing two median absolute deviations.
}
\usage{
plotBar(header, dataset, col4val, flag, plotTitle, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{col4val}{ a character string specifying the column whose intensity values will be used for the plot }
  \item{flag}{ 0, 1, or 2. 0 uses the data from the complete dataset, 1 generates one plot for each experiment, 2 generates one plot for each plate. }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
Saves the plots in pdf and png files named after the experiment name specified in the header concatenated with the \code{ plotTitle} and if applicable the experiment number and/or the plate number.


  When \code{flag} == 0, returns the plot name (\code{plotName}). 

  
  When \code{flag} == 1, returns a list containing:
  \item{plotName }{ The plot name }
  \item{minOfScreens}{ The number of the first experiment }
  \item{numOfScreens}{ The number of the last experiment }
  When \code{flag} == 2, returns a list containing: the plot name, a vector with the number of the first experiment and of the last experiment, and a vector with the number of the first plate and the number of the last plate.
}

\seealso{ \code{\link{ZScorePlot}}, \code{\link{ZScorePlotTwo}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

plotname <- plotBar(header, dataset, "SigIntensity", 0, "Data per well", 1)
}

\keyword{ device }