\name{LiWongRank}
\alias{LiWongRank}

\title{ Li Wong rank / invariant probeset normalization }
\description{
  Performs a Li Wong rank / invariant probeset normalization (see References).
}
\usage{
LiWongRank(header, dataset, listOfArgs)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{listOfArgs}{ a list containing:
  
  - a character string specifying the column whose values will be used for normalization

- a character string specifying the name of the dataset column to be used for the computation of the siRNA/gene ranks 
}
}
\details{
  For each plate type/layout in each experiment, generates a ranked list of siRNAs according to their intensity values. Only siRNAs occuring only once on the plate are allowed in the list. The normalization is performed only if all plate types have a maximum of 20 % of wells containing siRNAs occuring several times on the plate type. 

For each "unique" siRNA on a plate type, the variance of its ranks across plates is computed. A histogram of variances is plotted and allows the user to choose a threshold. A list of siRNAs with rank variances under the given threshold is then returned for each plate type so that the user can choose an siRNA to normalize the plate with.
}
\value{

Returns a list containing:
  \item{header }{the new header (with an added entry about the normalization procedure in the comments)}
  \item{dataset }{the new dataset with normalized values. The old values are saved in an extra column with the suffix ".old"}
}
\references{ C. Li and WH Wong. Model-based analysis of oligonucleotide arrays: model validation, design issues and standard error application. Genome Biol, 2(8):research0032.1-0032.11, 2001.

E. Schadt, C. Li, B. Ellis, and WH Wong. Feature Extraction and Normalization Algorithms for High-Density Oligonucleotide Gene Expression Array Data. J Cell Biochem Suppl, 37:120-125, 2001. }

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

normres <- LiWongRank(header, dataset, list("SigIntensity", "GeneName"))
newheader=normres[[1]]
newdataset=normres[[2]]
}

\keyword{ manip }
