\name{savepValVec}
\alias{savepValVec}

\title{ Save p-values to file }
\description{
  Saves a vector of p-values to a text file.
}
\usage{
savepValVec(pValVec, filename)
}

\arguments{
  \item{pValVec}{ a vector of p-values }
  \item{filename}{ the name of the text file to save the p-values to. }
}

\seealso{ \code{\link{Ttest}} }

\examples{
data(pValVec1, package="RNAither")

##for details on the generation of pValVec1, see the example of the Ttest function linked above.

savepValVec(pValVec1, "pvals_testfile1.txt")
}

\keyword{ datagen }
\keyword{ IO }