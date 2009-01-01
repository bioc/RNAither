\name{multTestAdjust}
\alias{multTestAdjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adjust p-values for multiple testing }
\description{
  Adjusts p-values for multiple testing.
}
\usage{
multTestAdjust(pValVec, adjustMethod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pValVec}{ a vector of p-values }
  \item{adjustMethod}{  one of the following: \code{ "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"}. For details type \code{?p.adjust}}
}

\value{
Returns a vector of corrected p-values. Can be integrated into a dataframe with the function \code{\link{incorporatepValVec}}.
}

\seealso{ \code{\link{incorporatepValVec}}, \code{\link{Ttest}} }

\examples{
data(pValVec1, package="RNAither")

##for details on the generation of pValVec1, see the example of the Ttest function linked above.

newpvalvec <- multTestAdjust(pValVec1,"fdr")
}
\keyword{ manip }
