\name{volcanoPlot}
\alias{volcanoPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Making a volcano plot }
\description{
  Makes a volcano plot of the data.
}
\usage{
volcanoPlot(header, dataset, col4plotx, col4ploty, col4anno, plotTitle, sigLevel, showPlot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{col4plotx}{ a character vector specifying the name of the column containing the intensity values, usually \code{ SigIntensity} }
  \item{col4ploty}{ a character vector specifying the name of the dataset column containing the corresponding p-values }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"}. }
  \item{plotTitle}{ the plot title }
  \item{sigLevel}{ the significance level for the p-value, indicating where a horizontal green line will be drawn }
  \item{showPlot}{ 0 or 1. 1 will open one or several plot windows in the R GUI, 0 will only save the plot(s) without opening windows }
}

\value{
Plots the intensity values against the negative decadic logarithm of the p-values. A green horizontal line is drawn at the specified significance level.

The plot is saved in a pdf and a png file named after the experiment name specified in the header concatenated with the \code{ plotTitle}.

The function returns the plot name.
}

\examples{
data(exampleHeader, package="RNAither")
data(pValVec1, package="RNAither")
data(scoredDataset1, package="RNAither")

##for details on the generation of pValVec1 and scoredDataset1, see \code{\link{Ttest}}

scoredHits1 <- hitselectionPval(scoredDataset1, pValVec1, "SigIntensity", "pValue.ttest_l", 0.05, 
"GeneName", "pvalue_testfile1.txt")

hitDataset1 <- scoredHits1[[1]]
hitvector1 <- scoredHits1[[2]]

volcano_name <- volcanoPlot(header, hitDataset1, "SigIntensity", "pValue.ttest_l", "GeneName", 
"Volcano Plot", 0.05, 1)
}

\keyword{ device }
\keyword{ univar }
