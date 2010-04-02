\name{numCellQualControl}
\alias{numCellQualControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quality control of the number of cells }
\description{
  Plots a histogram of the cell number per well and allows the user to set an upper and a lower threshold so as to exclude wells from the analysis.
}
\usage{
numCellQualControl(DataSetFile, nbLinesHeader, plotTitle)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataSetFile}{ a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{nbLinesHeader}{ typically 3 }
  \item{plotTitle}{ the plot title }
}

\value{
  Prints out the list of wells under and over the predefined thresholds in the shell.

  Saves a list of discarded siRNA values (if applicable) in a text file named after the experiment name specified in the header concatenated with either "numCellQualControl\_discarded\_higher.txt" or "numCellQualControl\_discarded\_lower.txt".

  Saves the histogram with the applied thresholds in a pdf file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

Overwrites the given \code{ DataSetFile} with the new dataset.
}

\seealso{ \code{\link{percCellQualControl}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")
saveDataset(header, dataset, "save_testfile1.txt")

numCellQualControl("save_testfile1.txt", 3, "Histogram of the number of cells")
}

\keyword{ device }
\keyword{ arith }
