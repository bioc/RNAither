\name{vennDiag}
\alias{vennDiag}

\title{ Plotting a Venn Diagram to compare hits }
\description{
  Plots a Venn Diagram of up to three binary hit vectors.
}
\usage{
vennDiag(header, listOfCols, listOfNames, plotTitle, showPlot)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{listOfCols}{ a list of binary hit vectors to compare }
  \item{listOfNames}{ a list of character strings for the annotation of the Venn Diagram }
  \item{plotTitle}{ the plot title }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
}

\value{
The plot is saved in a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns the plot name.
}

\examples{
data(exampleHeader, package="RNAither")

data(pValVec1, package="RNAither")
data(pValVec2, package="RNAither")
data(scoredDataset1, package="RNAither")
data(scoredDataset2, package="RNAither")

##for details on the generation of pValVec and scoredDataset, 
##see \code{\link{Ttest}} and \code{\link{MannWhitney}}

scoredHits1 <- hitselectionPval(scoredDataset1, pValVec1, "SigIntensity", "pValue.ttest_l", 0.05, 
"GeneName", "pvalue_testfile1.txt")

scoredHits2 <- hitselectionPval(scoredDataset2, pValVec2, "SigIntensity", "pValue.mannwhitney_l", 0.05, 
"GeneName", "pvalue_testfile2.txt")

hitvector1 <- scoredHits1[[2]]
hitvector2 <- scoredHits2[[2]]

plot_name <- vennDiag(header, list(hitvector1, hitvector2), list("t test", "Mann-Whitney test"), 
"Venn diagram", 1)
}

\keyword{ device }
