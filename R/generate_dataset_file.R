generateDatasetFile<-function(externalExperimentName, typeOfData, comments, 
outputFile, plateLayoutInternal, plateLayoutNCBI, nbRowsPerPlate, 
nbColsPerPlate, screenNb_pre, emptyWells, poorWells, controlCoordsOutput, 
backgroundValOutput, meanSignalOutput, SDmeanSignal, objNumOutput, 
cellNumOutput){

nbSpotsPerPlate<-nbRowsPerPlate*nbColsPerPlate

totalNumberOfSpots<-nbSpotsPerPlate*length(controlCoordsOutput)

Spotnumber<-rep(0, totalNumberOfSpots)
SpotType<-rep(0, totalNumberOfSpots)
Internal_GeneID<-rep("NA", totalNumberOfSpots)
GeneName<-rep("NA", totalNumberOfSpots)
SigIntensity<-rep(0, totalNumberOfSpots)
SDSIntensity<-rep(0, totalNumberOfSpots)
Background<-rep(0, totalNumberOfSpots)
LabtekNb<-rep(0, totalNumberOfSpots)
RowNb<-rep(0, totalNumberOfSpots)
ColNb<-rep(0, totalNumberOfSpots)
ScreenNb<-rep(0, totalNumberOfSpots)
NbCells<-rep(0, totalNumberOfSpots)
PercCells<-rep(0, totalNumberOfSpots)

    for (i in 1:length(controlCoordsOutput)){
    
        posCoords<-controlCoordsOutput[[i]][[1]]
        negCoords<-controlCoordsOutput[[i]][[2]]

        if (length(objNumOutput)>1){
            numObj<-objNumOutput[[i]]
        }else{
            if (sum(!is.na(objNumOutput))!=0){
                numObj<-objNumOutput[[i]]
            }else{
                numObj<-rep(as.double(NA_character_), nbSpotsPerPlate)
            }
        }

        if (length(cellNumOutput)>1){
            numCells<-cellNumOutput[[i]]
        }else{
            if (sum(!is.na(cellNumOutput))!=0){
                numCells<-cellNumOutput[[i]]
            }else{
                numCells<-rep(as.double(NA_character_), nbSpotsPerPlate)
            }
        }

        meanCytoSignal<-meanSignalOutput[[i]]

        SpotType_pre<-rep(2, nbSpotsPerPlate)
        
        if (!is.na(posCoords[1])){
            for (j in 1:length(posCoords)){
                SpotType_pre[posCoords[j]]<-1
            }
        }
        if (!is.na(negCoords[1])){
            for (k in 1:length(negCoords)){
                SpotType_pre[negCoords[k]]<-0
            }
        }
        if (!is.na(emptyWells[i])){
            SpotType_pre[emptyWells[i][[1]]]<-(-1)
        }
        if (!is.na(poorWells[i])){
            SpotType_pre[poorWells[i][[1]]]<-(-1)
        }
        
        ##compute vector for row nums:
        RowNb_pre<-rep(0, nbSpotsPerPlate)
        counter<-0
        for (l in seq(1, nbSpotsPerPlate, nbColsPerPlate)){
            counter<-counter+1
            RowNb_pre[l:(l+nbColsPerPlate-1)]<-rep(counter, nbColsPerPlate)
        }

        ##compute vector for col nums:
        ColNb_pre<-rep(0, nbSpotsPerPlate)
        counter<-0
        for (l in seq(1, nbSpotsPerPlate, nbColsPerPlate)){
            ColNb_pre[l:(l+nbColsPerPlate-1)]<-1:nbColsPerPlate
        }

        ##compute percent of obj recognized as cells:
        if (sum(is.na(numCells))<length(numCells) & sum(is.na(numObj))<length(numObj)){
            PercCells_pre<-as.vector(numCells)/as.vector(numObj)
        }else{
            PercCells_pre<-rep(as.double(NA_character_), nbSpotsPerPlate)
        }

        anfang<-1+(i-1)*nbSpotsPerPlate
        ende<-nbSpotsPerPlate*i
        Spotnumber[anfang:ende]<-1:nbSpotsPerPlate
        Internal_GeneID[anfang:ende]<-plateLayoutInternal[, i]
        GeneName[anfang:ende]<-plateLayoutNCBI[, i]
        SpotType[anfang:ende]<-SpotType_pre
        SigIntensity[anfang:ende]<-meanCytoSignal
        
        if (length(SDmeanSignal)>=i){
            if (length(SDmeanSignal[[i]]) == nbSpotsPerPlate){
                SDSIntensity[anfang:ende]<-SDmeanSignal[[i]]
            }else{
                SDSIntensity[anfang:ende]<-rep("NA", nbSpotsPerPlate)
            }
        }else{
            SDSIntensity[anfang:ende]<-rep("NA", nbSpotsPerPlate)
        }

        if (length(backgroundValOutput)>=i){
            if (length(backgroundValOutput[[i]]) == nbSpotsPerPlate){
                Background[anfang:ende]<-backgroundValOutput[[i]]
            }else{
                Background[anfang:ende]<-rep("NA", nbSpotsPerPlate)
            }
        }else{
            Background[anfang:ende]<-rep("NA", nbSpotsPerPlate)
        }
        LabtekNb[anfang:ende]<-rep(i, nbSpotsPerPlate)
        RowNb[anfang:ende]<-RowNb_pre
        ColNb[anfang:ende]<-ColNb_pre
        ScreenNb[anfang:ende]<-rep(screenNb_pre, nbSpotsPerPlate)
        NbCells[anfang:ende]<-numCells
        PercCells[anfang:ende]<-PercCells_pre

    }

    SpotType<-factor(SpotType, levels=c(-1, 0, 1, 2))

    dataFrame<-data.frame(Spotnumber, Internal_GeneID, GeneName, SpotType, 
    SigIntensity, SDSIntensity, Background, LabtekNb, RowNb, ColNb, ScreenNb, 
    NbCells, PercCells, stringsAsFactors=F)


    s1<-paste("external_experiment_name", externalExperimentName, sep=",")
    write.table(s1, file=outputFile, quote=F, col.names=F, row.names=F)
    
    s2<-paste("type_of_data", typeOfData, sep=",")
    write.table(s2, file=outputFile, quote=F, col.names=F, row.names=F, append=T)
    
    s3<-paste("comments", comments, sep=",")
    write.table(s3, file=outputFile, quote=F, col.names=F, row.names=F, append=T)

    write.table(dataFrame, file=outputFile, sep="\t", quote=F, append=T)
}


joinDatasets<-function(listOfDatasets){
    new<-rbind(listOfDatasets[[1]], listOfDatasets[[2]])
    if (length(listOfDatasets)>2){
        for (i in 3:length(listOfDatasets)){
            new<-rbind(new, listOfDatasets[[i]])
        }
    }
    rownames(new)<-1:nrow(new)
    invisible(new)
}

joinDatasetFiles<-function(listOfFiles, nbOfLinesInHeader, newHead, outputFile){
    listOfDatasets<-list("platzhalter")
    for (i in 1:length(listOfFiles)){
        data<-read.table(listOfFiles[[i]], skip=nbOfLinesInHeader)
        listOfDatasets[[i]]<-data
    }
    
    new<-rbind(listOfDatasets[[1]], listOfDatasets[[2]])
    if (length(listOfDatasets)>2){
        for (i in 3:length(listOfDatasets)){
            new<-rbind(new, listOfDatasets[[i]])
        }
    }
    rownames(new)<-1:nrow(new)
    
    write.table(newHead, file=outputFile, quote=F, col.names=F, row.names=F)
    write.table(new, file=outputFile, sep="\t", quote=F, append=T)

}
