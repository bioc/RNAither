\name{replicatesSpearmancor}
\alias{replicatesSpearmancor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the correlation coefficient betwenn replicates or experiments }
\description{
  Computes Spearman's rank correlation coefficient for each replicate - either inside each experiment, or between experiments.
}
\usage{
replicatesSpearmancor(header, dataset, flag, col4val, col4anno, fileNameSuffix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{flag}{ 1 or 2. 1 will compute the coefficient for a maximum of 3 replicates, for each experiment available in the dataset. 2 will summarize the replicates from each experiment with their root mean square and compute the correlation coefficient between experiments. }
  \item{col4val}{ a character string specifying the column whose values will be used to compute the correlation coefficient }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, e.g. \code{ "GeneName"} or \code{ "Internal_GeneID"} }
  \item{fileNameSuffix}{ a character string that will be used to name the output file containing a table with the correlation coefficients. }
}

\value{
For \code{ flag==1}, the correlation coefficients are printed out to the shell and saved in a text file named after the experiment name specified in the header concatenated with the character string filenamesuffix and "Spearmancor.txt".

For \code{ flag==2}, the correlation coefficients are printed out to the shell and saved in a text file named after the experiment name specified in the header concatenated with the character string filenamesuffix and "Spearmancor\_AllExp.txt".

The function returns a table containing the correlation coefficients.
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

replicatesSpearmancor(header, dataset, 1, "SigIntensity", "GeneName", "testfile1_")
}

\keyword{ arith }