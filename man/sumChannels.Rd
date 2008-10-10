\name{sumChannels}
\alias{sumChannels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summarize channels }
\description{
  Summarizes two channels, for example by computing their ratio.
}
\usage{
sumChannels(header, dataset, funName, colname4ch1, colname4ch2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{funName}{ the function used to summarize the two channels, for example \code{\link{divideChannels}} }
  \item{colname4ch1}{ a character string specifying the name of the dataset column containing the first channel }
  \item{colname4ch2}{ a character string specifying the name of the dataset column containing the second channel }
}

\details{
  The original dataset columns are saved as extra columns with the suffix \code{ ".old"} by the function \code{ \link{saveOldIntensityColumns}}.
}

\value{ A list containing:

  \item{header }{the header with an entry about the channel summarization added in the comments section}
  \item{newDataset }{the new dataset}
}

\seealso{ \code{\link{eraseDataSetColumn}}, \code{\link{divideChannels}}, \code{\link{saveOldIntensityColumns}} }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

newdataset=sumChannels(header, dataset, divideChannels, "SigIntensity", "NbCells")
}

\keyword{ manip }
