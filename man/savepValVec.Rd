\name{savepValVec}
\alias{savepValVec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Save p-values to file }
\description{
  Saves a vector of p-values to a text file.
}
\usage{
savepValVec(pValVec, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pValVec}{ a vector of p-values }
  \item{filename}{ the name of the text file to save the p-values to. }
}

\examples{
data(pValVec1, package="RNAither")

##for details on the generation of pValVec1, see \code{\link{Ttest}}

savepValVec(pValVec1, "pvals_testfile1.txt")
}

\keyword{ datagen }
\keyword{ IO }