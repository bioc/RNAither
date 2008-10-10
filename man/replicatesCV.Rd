\name{replicatesCV}
\alias{replicatesCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the correlation of variation (CV) }
\description{
  Computes the correlation of variation as defined in Tseng et al. (see References)
}
\usage{
replicatesCV(header, dataset, PlotTitle, col4val, col4anno, plotDesign, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{PlotTitle}{ the plot title }
  \item{col4val}{ a character string specifying the column whose values will be used to compute the correlation of variation }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"} }
  \item{plotDesign}{ 1 or 2. 1 will generate one window containing all plots, 2 will generate a series of plots. }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
The correlation of variation of an siRNA is defined as the standard deviation of its values divided by their mean.

The function generates a plot of the average intensity against the CV for each experiment. The plot will be saved as a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ PlotTitle}.

The function returns a list containing:

  \item{histoName }{the plotname}
  \item{minOfScreens }{the number of the first experiment}
  \item{numOfScreens }{the number of the last experiment}
}
\references{ G. C. Tseng et al. Issues in cDNA microarray analysis: quality filtering, channel normalization, models of variations and assessment of gene effects. Nucleic Acids Res, 29(12): 2549-2557, 2001. }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

replicatesCV(header, dataset, "Correlation of Variation versus Mean Intensity", 
"SigIntensity", "GeneName", 1, 0)
}

\keyword{ arith }
