\name{findReplicates}
\alias{findReplicates}

\title{ Find all replicates of a certain siRNA/gene in a dataset }
\description{
  Returns which lines in the dataset correspond to a given siRNA/gene ID.
}
\usage{
findReplicates(dataset, whichCol, replicateID)
}

\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{whichCol}{ a character string specifying the name of the column containing the ID, either \code{ Internal_GeneID} or \code{ GeneName} }
  \item{replicateID}{ the siRNA/gene ID of interest }
}

\value{
An integer vector containing the indexes in the main dataset of all wells corresponding to a given siRNA/gene ID
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

Indexes <- findReplicates(dataset, "GeneName", "CPSF1")
}

\keyword{ methods }
