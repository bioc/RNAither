\name{ZScorePlot}
\alias{ZScorePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot normalized intensity values per well }
\description{
  Plots the normalized intensity values for each well, together with a black line showing the mean, two green lines showing the standard deviation, and two red lines showing 2 standard deviations.
}
\usage{
ZScorePlot(header, dataset, flag, col4plot, col4anno, plotTitle, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{flag}{ either 1 or 2. 1 if the dataset contains values per well, 2 if the dataset contains summarized values for each siRNA (e.g. a dataset summarized with \code{\link{summarizeReps}}). }
  \item{col4plot}{ a character string specifying the column whose values will be used for the plot }
  \item{col4anno}{ a character string specifying the column that will be used for the plot annotation }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{  0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
Plots the normalized intensity values for each well, together with a black line showing the mean, and two red lines showing 2 standard deviations. Clicking on the points shows the gene/siRNA name.

The plot is saved as a pdf and a png file named after the experiment name specified in the header concatenated with the \code{plotTitle}.

The function returns the plot name.
}

\seealso{ \code{\link{plotBar}}, \code{\link{ZScorePlotTwo}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normedvals <- ZScore(header, dataset, list("SigIntensity", 1))
ZScorePlot(normedvals[[1]], normedvals[[2]], 1, "SigIntensity", "GeneName",
"Normed intensity values per well", 1)
}

\keyword{ device }
