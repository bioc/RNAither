\name{joinDatasetFiles}
\alias{joinDatasetFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Join dataset files }
\description{
  Merges two or more dataset files into one, with one common header.
}
\usage{
joinDatasetFiles(listOfFiles, nbOfLinesInHeader, newHead, outputFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listOfFiles}{ a list of the names of the files to join }
  \item{nbOfLinesInHeader}{ typically 3 }
  \item{newHead}{ the new header }
  \item{outputFile}{ the name of the file to save the header and concatenated dataset in }
}

\seealso{ \code{\link{generateDatasetFile}}, \code{\link{joinDatasets}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")
saveDataset(header, dataset, "save_testfile1.txt")

header[[1]] <- "external_experiment_name,Test screen"
header[[2]] <- "comments,contains twice Screen Nb 1"

joinDatasetFiles(list( "save_testfile1.txt", "save_testfile1.txt"), 3, header, 
"concatenated_testfile.txt")
}

\keyword{ datagen }
