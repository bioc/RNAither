\name{hitselectionZscore}
\alias{hitselectionZscore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Selecting hits according to ZScores }
\description{
  Selects significant genes according to their ZScore.
}
\usage{
hitselectionZscore(dataset, col4zscore, col4sel, thresh, flag, flag2, col4anno, 
sumFunc, file4hits)
}

\arguments{
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  \item{col4zscore}{ a character vector specifying the name of the column containing the ZScores, usually \code{SigIntensity} }
  \item{col4sel}{ a character vector specifying the name of the new dataset column where hits will be stored }
  \item{thresh}{ the threshold for the ZScores. The interpretation depends on the choice of the parameter \code{ flag2}. }
  \item{flag}{ 1 or 2. 1 means the ZScores are kept per well, 2 that they are summarized according to the parameter \code{ sumFunc}. }
  \item{flag2}{ 1, 2 or -2. 
  
  If 1 is chosen and \code{ thresh} == n, then the n greatest Zscores are chosen as hits. 
  
  If 1 is chosen and \code{ thresh} == -n, then the n smallest Zscores are chosen. 
  
  If 1 is chosen and \code{ thresh} == 0, all ZScores are chosen and written to the output file. 
  
  If 2 is chosen, all Zscores greater than or equal to \code{ thresh} are chosen. 
  
  If -2 is chosen, all Zscores smaller than or equal to \code{ thresh} are chosen.  }
  \item{col4anno}{ a character string specifying the name of the dataset column to be used to define the replicate, for example \code{ "GeneName"} or \code{ "Internal_GeneID"} }
  \item{sumFunc}{ the function used to summarize ZScore values, e.g. \code{ mean} or \code{ median}. }
  \item{file4hits}{ the name of the file to store the results in }
}
\details{
  If \code{ flag2} == -2, and there are no ZScores under the defined threshold \code{ thresh}, the threshold is increased to \code{ min(ZScores)}. 
  
  If \code{ flag2} == 2, and there are no ZScores over the defined threshold \code{ thresh}, the threshold is increased to \code{ max(ZScores)}.
}
\value{

A list containing: 

  \item{dataset }{the dataset with an added column defining the hits in the form of a binary vector}
  \item{hitVector }{the binary vector itself}
  \item{thresh }{the threshold for the ZScores}  
  ZScores are stored in a text output file.
}
\references{ N. Malo et al. Statistical practice in high-throughput screening data analysis. Nature Biotech, 24(2): 167-175, 2006. }

\seealso{ \code{\link{hitselectionPval}}, \code{\link{hitselectionZscorePval}} }
\examples{
data(scoredDataset1, package="RNAither")
data(pValVec1, package="RNAither")

##for details on the generation of pValVec1 and scoredDataset1, see \code{\link{Ttest}}

scoredHits1 <- hitselectionZscore(scoredDataset1, "SigIntensity", "Zscore_hits", -10, 
2, 1, "GeneName", median, "Zscores_testfile1.txt")

newdataset <- scoredHits1[[1]]
hitvector <- scoredHits1[[2]]
}
\keyword{ arith }