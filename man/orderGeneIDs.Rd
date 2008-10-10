\name{orderGeneIDs}
\alias{orderGeneIDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Order a dataset }
\description{
  Orders dataset according to one of its columns.
}
\usage{
orderGeneIDs(dataset, ID1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{ID1}{ a character string specifying the name of the column according to which the dataset will be sorted }
}

\value{
An R data frame ('dataset') ordered according to its values in the specified column.
}

\seealso{ \code{\link{order}} }
\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

orderedDataset=orderGeneIDs(dataset,"SigIntensity")
}

\keyword{ arith }
