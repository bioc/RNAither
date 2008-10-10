\name{mainAnalysis}
\alias{mainAnalysis}

\title{ Wrapper function for full automated analysis }
\description{
  Performs a standard analysis of the data (quality and statistics) from a dataset file.
}
\usage{
mainAnalysis(header, dataset, flagForSameExp, listOfNormalizations, listOfArgs4norm, 
listOfStatTests, listOfArgs4stat, multTestAdj, hitScoringVec1, hitScoringVec2, 
posNegFlag, flag4Gsea, vecOfChannels, whichOnto)
}

\arguments{
  \item{header}{ the header of a dataset file generated with \code{\link{generateDatasetFile}} }
  
  \item{dataset}{ an R data frame generated with \code{\link{generateDatasetFile}} }
  
  \item{flagForSameExp}{ either 0 or 1. If 1, all experiments defined in the column \code{ ScreenNb} in the dataset file must have the same design (same type and same number of replicates - exact plate layout is irrelevant) so that Spearman's correlation coefficient can be computed between experiments (each with summarized replicates) }
  
  \item{listOfNormalizations}{ a list of the normalization function to apply. Can be \code{\link{LiWongRank}}, \code{\link{varAdjust}}, \code{\link{divNorm}}, \code{\link{quantileNormalization}}, \code{\link{BScore}}, \code{\link{ZScore}}, \code{\link{ZScorePerScreen}}, \code{\link{subtractBackground}}, \code{\link{lowessNorm}}, \code{\link{controlNorm}} }
  
  \item{listOfArgs4norm}{ a list containing, for each element of \code{ listofnormalizations}, the arguments to be passed on }
  
  \item{listOfStatTests}{ a list of the statistical tests to perform. Can be \code{\link{Ttest}}, \code{\link{MannWhitney}}, \code{\link{RankProduct}} }
  
  \item{listOfArgs4stat}{ a list containing, for each element of \code{ listofstattests}, the arguments to be passed on }
  
  \item{multTestAdj}{ indicates the p-value correction for multiple testing - one of \code{ "holm"}, \code{ "hochberg"}, \code{ "hommel"}, \code{ "bonferroni"}, \code{ "BH"}, \code{ "BY"}, \code{ "fdr"}, or \code{ "none"} (Type \code{?p.adjust} for details))}
  
  \item{hitScoringVec1}{ a vector of length 3 indicating (in that order): 
  
  - scoring according to p-value (0: no, 1: yes)
  
  - scoring according to ZScore with ZScore < threshold (0: no, 1: yes), or according to ZScore < threshold and p-value < hitScoringVec2[1] (2)
  
  - scoring according to ZScore with ZScore > threshold (0: no, 1: yes), or according to ZScore > threshold and p-value < hitScoringVec2[1] (2). 
  
  If hitScoringVec1[2] or hitScoringVec1[3] are equal to 2, hitScoringVec1[1] must be equal to one, otherwise p-values will not be computed. }
  
  \item{hitScoringVec2}{ a vector of length 3 indicating the thresholds for hitscoringvec1 }
  
  \item{posNegFlag}{ either 0 (no controls available) or 1 (controls available) }
  
  \item{flag4Gsea}{ Can be:
  
  - either 0: No GSEA analysis is performed
  
  - or 1: A GSEA analysis is performed for each hit scoring method
  
  - or a binary vector that allows to choose which hit scoring method(s) will be used for a GSEA analysis. Hit scoring methods are sorted as follows: first, hits are scored according to the p-values of each test specified in \code{ listOfStatTests }. Then, if the option of scoring hits according to p-values and Intensities is chosen (see \code{ hitScoringVec1 }, for each test, a hit vector is generated. Finally, if the option of scoring hits according to Intensities only is chosen, hit vectors are generated for this option.
  }
  
  \item{vecOfChannels}{ a character vector containing the names of the channels to be used for quality plots, for example \code{ "SigIntensity"} or \code{ "NbCells"} }
  
  \item{whichOnto}{ one of the three GO hierarchies: \code{ "biological_process" }, \code{ "molecular_function" } or \code{ "cellular_component" } - used for the GSEA analysis }
}

\value{
Generates the html output files \code{ index.html} and \code{ indexnorm.html} containing the quality analysis of raw and normalized data, respectively, and \code{ stats.html}, containing the statistical analysis. If several normalization methods are applied, an \code{indexnorm} file is generated after each.
}

\examples{
data(exampleHeader, package="RNAither")
data(exampleDataset, package="RNAither")

mainAnalysis(header, dataset, 0, list(controlNorm), list(list(1, 0, "SigIntensity", 1)), 
list(Ttest, MannWhitney), list(list("l", 1, "SigIntensity", "GeneName"), 
list("l", 1, "SigIntensity", "GeneName")), "none", c(1, 0, 0), c(0.05, 0, 0), 1, 
1, c("SigIntensity", "NbCells"), "biological_process")

}

\keyword{ dplot }
\keyword{ device }
\keyword{ arith }
\keyword{ manip }
\keyword{ database }
\keyword{ methods }
\keyword{ univar }
\keyword{ datagen }
\keyword{ misc }
\keyword{ utilities }