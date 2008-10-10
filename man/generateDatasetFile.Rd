\name{generateDatasetFile}
\alias{generateDatasetFile}

\title{ Generate Dataset File }
\description{
  Generates a text file containing all experimental data. Needed for all subsequent analysis functions.
}
\usage{
generateDatasetFile(externalExperimentName, typeOfData, comments, outputFile, 
plateLayoutInternal, plateLayoutNCBI, nbRowsPerPlate, nbColsPerPlate, screenNb_pre, 
emptyWells, poorWells, controlCoordsOutput, backgroundValOutput, meanSignalOutput, 
SDmeanSignal, objNumOutput, cellNumOutput)
}

\arguments{
  \item{externalExperimentName}{ A character string specifying the experiment name, e.g. "Johns Experiment Nb. 1" }
  
  \item{typeOfData}{  A character string specifying the type of data, e.g. "364 well plate data for virus screens" }
  
  \item{comments}{ A character string specifying comments. NA if not available. }
  
  \item{outputFile}{ A character string specifying the name of the text file containing the dataset. }
  
  \item{plateLayoutInternal}{ A matrix of internal siRNA IDs specifying their position on the plate (row-wise). Each column of the matrix stands for one plate. }
  
  \item{plateLayoutNCBI}{ A matrix of gene names specifying their position on the plate (row-wise). Each column of the matrix stands for one plate. }
  
  \item{nbRowsPerPlate}{ The number of rows per plate }
  
  \item{nbColsPerPlate}{ The number of columns per plate }
  
  \item{screenNb_pre}{ The screen/experiment number }
  
  \item{emptyWells}{ A list containing, for each plate, an integer vector of the positions of empty wells. NA if there are no empty wells on the plate. }
  
  \item{poorWells}{ A list containing, for each plate, an integer vector of the positions of wells that, for a certain reason, should not be taken into account during the analysis. NA if there are no such wells on the plate. }
  
  \item{controlCoordsOutput}{ A list containing, for each plate, a list of integer vectors specifying the positions of positive (first element in sublist) and negative (second element in sublist) controls. NA if there are no positive/negative controls on the plate. }
  
  \item{backgroundValOutput}{ A list containing, for each plate, a vector of background values per well }
  
  \item{meanSignalOutput}{ A list containing, for each plate, a vector of intensity values for each well }
  
  \item{SDmeanSignal}{ A list containing, for each plate, a vector of standard deviations of intensity values for each well }
  
  \item{objNumOutput}{ A list containing, for each plate, a vector of the number of identified objects for each well }
  
  \item{cellNumOutput}{ A list containing, for each plate, a vector of intensity values for each well, e.g. a vector of the number of identified cells for each well. }
}
\details{
  Positions on plates are specified with one integer only. For example, the position of the well in row 2 and column 5 is \code{ (RowNo-1)*(Number of columns on plate)+ColNo}.
}
\value{
  The function generates a text file consisting of a header and a 'dataset'. The header contains the experiment description (\code{ExternalExperimentName}, \code{TypeOfData} and \code{Comments}). The dataset is an R data frame, each row corresponding to one well, with the following columns:
  
  \item{Spotnumber}{ The position of the well on the plate}
  
  \item{Internal_GeneID}{ The ID of the siRNA}
  
  \item{GeneName}{ The gene name}
  
  \item{SpotType}{ Can be -1, 0, 1 or 2.
  
  Type -1 wells (e.g. emtpy wells, wells with poor quality) are not considered in subsequent analyses but are kept in the data set for the sake of completeness. 
  
  Type 0 wells correspond to negative controls, type 1 wells to positive controls. 
  
  Type 2 wells correspond to the standard data wells.}
  
  \item{SigIntensity}{ The signal intensity (channel 1)}
  
  \item{SDSIntensity}{ The standard deviation of the signal intensity, if available}
  
  \item{Background}{ The background per well, if available}
  
  \item{LabtekNb}{ The plate number}
  
  \item{RowNb}{ The row number}
  
  \item{ColNb}{ The column number}
  
  \item{ScreenNb}{ The screen number}
  
  \item{NbCells}{ E.g. the number of cells identified in the well (channel 2)}
  
  \item{PercCells}{ The ratio (number of identified cells)/(number of identified objects)}
}
\seealso{ \code{\link{joinDatasetFiles}}, \code{\link{joinDatasets}} }

\examples{

##gene names
plateLayout1 <- c("test1", "empty", "test3", "test4", "test5", 
"test6", "test7", "empty", "test9", "test10", "test11", "test12")

plateLayout2 <- c("test1", "test2", "test3", "test4", "test5", 
"test6", "test7", "test8", "test9", "test10", "test11", "test12")

plateLayout <- cbind(plateLayout1, plateLayout2)

emptyWells <- list(c(2, 8), NA_integer_)
##the first plate has two empty wells at position 2 and 8,
##the second plate does not have any empty wells

poorWells <- NA_integer_
##no wells of poor quality

controlCoordsOutput <- list(list(NA_integer_, NA_integer_), list(NA_integer_, c(9,10)))
##the first plate does not have any control siRNAs,
##the second plate has two negative controls at position 9 and 10

backgroundValOutput<-NA_integer_
##no background signal intensities available

sigPlate1<-c(2578, NA_integer_, 3784, 3784, 2578, 5555, 5555, NA_integer_, 8154, 2578, 3784, 2578)
sigPlate2<-c(8154, 3784, 5555, 3784, 11969, 2578, 1196, 5555, 17568, 2578, 5555, 2578)
##the signal intensities on the plates

meanSignalOutput<-list(sigPlate1, sigPlate2)

SDmeansignal<-NA_integer_
##no standard deviation available

objnumOutput<-NA_integer_
##no cell count available

cellnumOutput<-NA_integer_

generateDatasetFile("First test screen", "RNAi in virus-infected cells", 
NA_character_, "testscreen_output.txt", plateLayout, plateLayout, 3, 4, 
1, emptyWells, poorWells, controlCoordsOutput, backgroundValOutput, 
meanSignalOutput, SDmeansignal, objnumOutput, cellnumOutput)

##load the dataset into R:
header<-readLines("testscreen_output.txt",3)
dataset<-read.table("testscreen_output.txt", skip=3, colClasses=c(NA, NA, NA, NA, 
"factor", NA, NA, NA, NA, NA, NA, NA, NA, NA), stringsAsFactors=FALSE)
}
\keyword{ datagen }
