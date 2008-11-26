\name{DRQualControl}
\alias{DRQualControl}

\title{ Computing the dynamic range }
\description{
  Computes the dynamic range per plate for a complete dataset file and plots the results.
}
\usage{
DRQualControl(header, data, nbLinesHeader, channel, plotTitle, showPlot)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{data}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{nbLinesHeader}{ typically 3 }
  \item{channel}{ A character string specifying the name of the column containing the values for computing the dynamic range, e.g. \code{ "SigIntensity"} }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
}

\value{
Returns the dynamic range for each plate in the shell and saves them in a text file. The name of the text file will be the concatenation of the experiment name specified in the header and the character string "DR.txt".

Shows a plot of the dynamic range values and saves it as a pdf file under the experiment name specified in the header concatenated with the function argument \code{ plotTitle}.
}

\references{ M. Boutros, L. Bras, and W. Huber. Analysis of cell-based RNAi screens. \textit{Genome Biol}, 7(7): R66, 2006. }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

DRQualControl(header, dataset, 3, "SigIntensity", "DR per plate", 1)
}

\keyword{ device }

