\name{SNRQualControl}
\alias{SNRQualControl}

\title{ Computing the SNR }
\description{
  Computes the signal to noise ratio for all data, per experiment and per plate for a complete dataset file and plots histograms of the results.
}
\usage{
SNRQualControl(dataSetFile, nbLinesHeader, channel, noise, plotTitle, showPlot)
}

\arguments{
  \item{dataSetFile}{ a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{nbLinesHeader}{ typically 3 }
  \item{channel}{ a character string specifying the name of the column containing the values for computing the SNR, e.g. \code{"SigIntensity"} }
  \item{noise}{ A character string specifying the name of the column containing the values for computing the SNR, e.g. \code{"Background"} }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
}

\value{
Shows histogram plots of the SNR for the whole dataset file, per experiment and per plate and saves them in a pdf file. The name of the file will be the concatenation of the experiment name specified in the header and the function argument \code{ plotTitle}.
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")
saveDataset(header, dataset, "save_testfile1.txt")

SNRQualControl("save_testfile1.txt", 3, "SigIntensity", "Background", "SNR", 1)
}

\keyword{ arith }

