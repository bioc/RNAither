\name{summarizeReps}
\alias{summarizeReps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate a new dataset with summarized replicates }
\description{
  Generates a new dataset with summarized replicates.
}
\usage{
summarizeReps(data, funSum, col4val, col4anno, cols2del)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{funSum}{ a function used to summarize the values of a replicate, e.g. \code{ mean}, \code{ median}, \code{\link{rms}}, \code{\link{trim}}, \code{ max}, \code{ min}, \code{\link{closestToZero}}, \code{\link{furthestFromZero}}, ... }
  \item{col4val}{ a character vector (containing for example \code{ "SigIntensity"}, \code{ Background}, \code{ NbCells}, \code{ PercCells}, ...) specifying the columns that will be summarized by \code{ funSum} }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"} }
  \item{cols2del}{ a character vector containing the columns to delete, for example \code{ "SDSIntensity"} }
}
\details{
  All columns containing replicate values will be summarized by \code{ funSum}. For all columns containing positions, screen numbers, plate numbers, etc., all information for different replicates will be kept, comma-separated. All columns containing standard deviations of channels should be specified in \code{ colnames2delete}.
}
\value{
  Returns the summarized dataset.
}

\seealso{ \code{\link{summarizeRepsNoFiltering}}, \code{\link{eraseDataSetColumn}}, \code{\link{generateReplicateMat}}, \code{\link{generateRepMatNoFilter}}, \code{mean}, \code{median}, \code{\link{rms}}, \code{\link{trim}}, \code{max}, \code{min}, \code{\link{closestToZero}}, \code{\link{furthestFromZero}} }

\examples{
data(exampleDataset, package="RNAither")

colname4val <- c("SigIntensity", "Background", "NbCells", "PercCells")
summarizeddataset <- summarizeReps(dataset, mean, colname4val, "GeneName", "SDSIntensity")
}

\keyword{ arith }