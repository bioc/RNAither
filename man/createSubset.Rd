\name{createSubset}
\alias{createSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creating a subset of a dataset according to a certain column value }
\description{
  Creates a subset of a dataset containing all wells/lines having a certain value in a specified column.
}
\usage{
createSubset(dataset, listIDs, equalTo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listIDs}{ a character string and one of the following: \code{ Spotnumber, Internal_GeneID, GeneName, SpotType, SigIntensity, SDSIntensity, Background, LabtekNb, RowNb, ColNb, ScreenNb, NbCells, PercCells, ...} }
  \item{equalTo}{ A value or character string specifying the value in the chosen column, e.g. all wells on plate 2 }
}
\value{
A subset of the dataset containing only the wells/lines having a certain value in a specified column.
}

\seealso{ \code{\link{indexSubset}} }
\examples{
data(exampleDataset, package="RNAither")

subset <- createSubset(dataset, dataset$LabtekNb, 2)
}
\keyword{ arith }
