\name{ZPRIMEQualControl}
\alias{ZPRIMEQualControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computing the Z' factor }
\description{
  Computes the Z' factor per plate for a complete dataset file and plots the results.
}
\usage{
ZPRIMEQualControl(header, data, channel, plotTitle, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{data}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{channel}{ a character string specifying the name of the column containing the values for computing the Z' factor, e.g. \code{"SigIntensity"} }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows. }
}

\value{
Returns the Z' values in the shell for each plate and saves them in a text file. The name of the text file will be the concatenation of the experiment name specified in the header and the character string "Z'Scores.txt".

Shows a plot of the Z' factor values and saves it as a png and a pdf file under the experiment name specified in the header concatenated with the function argument \code{ plotTitle}.


The function returns a list containing:

  \item{plotName }{the plot name}
  \item{ZPrimeTabelle }{table containing the Z' values}
}
\references{ J. Zhang, T. Chung, and K. Oldenburg. A simple statistical parameter for use in evaluation and validation of high throughput screening assays. J Biomol Screen, 4:67-73, 1999. }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

res <- ZPRIMEQualControl(header, dataset, "SigIntensity", "Z' factors per plate", 1)
zprime_plot <- res[[1]]
zprime_table <- res[[2]]
}

\keyword{ arith }