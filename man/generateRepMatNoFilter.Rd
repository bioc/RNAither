\name{generateRepMatNoFilter}
\alias{generateRepMatNoFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate a matrix of replicates (II) }
\description{
  Generates a matrix out of a dataset, each row corresponding to an siRNA/gene ID, each column to a channel value or its index in the dataset.
}
\usage{
generateRepMatNoFilter(data, minNbReps, IndexOrInt, col4val, col4anno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{minNbReps}{ set to 2 if you want to exclude replicates occurring only once in the dataset, otherwise 1. }
  \item{IndexOrInt}{ a character string - either \code{ "Index"} or \code{ "Intensities"} - specifying which values are to be contained in the output matrix. }
  \item{col4val}{ a character string specifying the name of the dataset column to be used for the values of the output matrix (if \code{ IndexOrIntensities} is set to \code{ "Intensities"}), for example \code{ "SigIntensity"} or \code{ "NbCells"} }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used for the output matrix' rows, for example \code{ "GeneName"} or \code{ "Internal_GeneID"}. }
}
\details{
  The function will not omit values or indexes of lines/wells with spot type -1. If you want to omit those, use \code{ generateReplicatematrix}.
}
\value{
A matrix with each row corresponding to an siRNA/gene ID (as reflected in rownames), each column to a channel value or its index in the dataset. Missing values (in case of different number of replicates occuring for different siRNAs/genes) are set to NA.
}

\seealso{ \code{\link{generateReplicateMat}} }
\examples{
data(exampleDataset, package="RNAither")

replicatematrix <- generateRepMatNoFilter(dataset, 2, "Index", "SigIntensity", "GeneName")
}

\keyword{ arith }
