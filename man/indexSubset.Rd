\name{indexSubset}
\alias{indexSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Saving the indexes of a subset in the main dataset }
\description{
  Used together with \code{\link{createSubset}}, returns the indexes in the main dataset of the wells chosen as a subset by the previous call of \code{\link{createSubset}}.
}
\usage{
indexSubset(listIDs, equalTo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listIDs}{ a character string and one of the following: \code{ Spotnumber, Internal_GeneID, GeneName, SpotType, SigIntensity, SDSIntensity, Background, LabtekNb, RowNb, ColNb, ScreenNb, NbCells, PercCells, ...} }
  \item{equalTo}{ A value or character string specifying the value in the chosen column, e.g. all wells on plate 2 }
}

\value{
An integer vector containing the indexes in the main dataset of the wells chosen as a subset by the previous call of \code{\link{createSubset}}.
}

\seealso{ \code{\link{createSubset}} }
\examples{
data(exampleDataset, package="RNAither")

subset <- createSubset(dataset, dataset$LabtekNb, 2)
indexOfSubsetInDataset <- indexSubset(dataset$LabtekNb, 2)
}
\keyword{ arith }
